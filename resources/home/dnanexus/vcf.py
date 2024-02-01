###
# functions to apply to vcf from eggd_vep
###
import os
import subprocess

import pandas as pd
from collections import defaultdict
from pathlib import Path
from pysam import VariantFile

# Get path of parent directory
ROOT_DIR = Path(__file__).absolute().parents[1]

def decompress(input_file):
    """
    Parameters
    -----------
    vcf : file
    path to file to be decompressed

    returns
    --------
    decompressed vcf file

    """

    print(f"Calling bgzip -d on {file}")

    output_file = subprocess.run(['gzip', '-d', file])

    return output_file


def bcftools_pre_process(input_vcf) -> str:
    """
    Decompose multiple transcript annotations to individual records, and split
    VEP CSQ string fields to individual INFO keys. Adds a 'CSQ_' prefix to
    these fields to stop potential conflicts with existing INFO fields. Then
    strips the original INFO/CSQ fields

    Parameters
    ----------
    input_vcf : file
        path to VCF file to be split

    Returns
    -------
    output_vcf : str
        name of VCF split by bcftools

    Raises
    ------
    AssertionError
        Raised when non-zero exit code returned by bcftools or output VCF
        is truncated (fewer variants than were in input VCF)
    """

    print(
        f"Splitting necessary fields from {input_vcf} "
        "with bcftools +split-vep"
    )
    output_vcf = f"{Path(input_vcf).stem.split('.')[0]}.split.vcf"

    # Check total rows before splitting out columns
    pre_split_total = subprocess.run(
        f"zgrep -v '^#' {input_vcf} | wc -l",
        shell=True,
        capture_output=True
    )
    # Split out all fields from the CSQ string and name them with
    # 'CSQ_{field}' as separate INFO fields. Remove original CSQ string
    cmd = (
        f"bcftools +split-vep --columns - -a CSQ -Ou -p 'CSQ_'"
        f" -d {input_vcf} | bcftools annotate -x INFO/CSQ -o {output_vcf}"
    )

    output = subprocess.run(cmd, shell=True, capture_output=True)

    assert output.returncode == 0, (
        f"\n\tError in splitting VCF with bcftools +split-vep. VCF: {input_vcf}"
        f"\n\tExitcode:{output.returncode}"
        f"\n\t{output.stderr.decode()}"
    )

    # Check total rows after splitting
    post_split_total = subprocess.run(
        f"zgrep -v '^#' {output_vcf} | wc -l",
        shell=True,
        capture_output=True
    )

    pre_split = pre_split_total.stdout.decode()
    post_split = post_split_total.stdout.decode()

    print(
        f"Total lines before splitting: {pre_split}"
        f"Total lines after splitting: {post_split}"
    )

    assert post_split >= pre_split, (
        "Count of variants following bcftools +split-vep is fewer than the "
        "input VCF"
    )

    return output_vcf


def read_in_vcf(vcf_file):
    """
    Read in the VCF file with the pyvcf package and add a new header line
    for MANE field

    Parameters
    ----------
    vcf_file : string
        name of the annotated VCF file 

    Returns
    -------
    vcf_contents: pyvcf.Reader object
        contents of the VCF file as a pyvcf object
    sample_name : str
        the full name of the sample
    """
    print(f"Reading in the split VCF {vcf_file} with pySAM")

     # Read in and create pysam object of the VCF
    vcf_contents = VariantFile(vcf_file, 'r')

    # Get the name of the sample from the VCF
    sample_name = list(vcf_contents.header.samples)[0]

    # Add MANE as INFO field
    vcf_contents.header.info.add(
        "MANE", "0", "Flag",
        "Check if Matched Annotation from NCBI and EMBL-EBI (MANE)")

    return vcf_contents, sample_name


def add_MANE_field(vcf_contents, transc) -> dict:
    """
    Add MANE INFO field to each variant 

    Parameters
    ----------
    vcf_contents : pysam.VariantFile object
        pysam object containing all the VCF's info
    panel_dict : dict
        default dict with gene symbol as key and gene info as val

    Returns
    -------
    gene_variant_dict : dict
        dictionary of each gene with value as a list of all variants in that
        gene (plus additional INFO field)
    """
    # Add each variant in a gene/entity to a dict, with gene as key and list
    # of variants as value
    gene_variant_dict = defaultdict(list)
    for variant_record in vcf_contents:
        gene = variant_record.info['CSQ_Feature'][0]
        gene_variant_dict[gene].append(variant_record)

    # For each gene, check whether certain gene info is available
    for gene, variant_list in gene_variant_dict.items():
        gene_present = gene_mane = False
        if gene in transc:
            gene_present = True
            gene_mane = "MANE"

        # Iterate over all of the variants called in that gene
        # If gene present in mane_list and gene_mane is present,
        # add the MANE we've taken from PanelApp to variant info
        # otherwise set it to unknown
        for variant in variant_list:
            if all([gene_present, gene_mane]):
                variant.info['MANE'] = 1
            else:
                variant.info['MANE'] = 0

    return gene_variant_dict


def write_out_flagged_vcf(flagged_vcf, gene_variant_dict, vcf_contents):
    """
    Write out each variant record to VCF using pysam

    Parameters
    ----------
    flagged_vcf : str
        Name of the VCF to be written out with flags added
    gene_variant_dict : dict
        dictionary with each gene and value as all of the variants in that gene
        as list
    vcf_contents : pysam.VariantFile object
        the contents of the VCF as a pysam object
    """
    print(f"Writing out flagged variants to VCF: {flagged_vcf}")

    with VariantFile(flagged_vcf, 'w', header=vcf_contents.header) as out_vcf:
        # For each gene, write out each variant to VCF with extra INFO field
        for gene, variant_list in gene_variant_dict.items():
            for variant in variant_list:
                out_vcf.write(variant)


def check_written_out_vcf(
        original_vcf_contents, gene_variant_dict, flagged_vcf
    ):
    """
    Check that the VCF file written out is exactly the same as the header
    and dict of pysam variants that was meant to be written

    Parameters
    ----------
    original_vcf_contents : pysam.VariantFile object
        the pysam object which was to be written to file
    gene_variant_dict : dict
        dictionary with each gene and value as all of the variants in that gene
        as list
    flagged_vcf : file
        the VCF that was written out to be read back in with pysam
    Raises
    ------
    AssertionError
        Raised if the contents of the VCF which was to be written out does
        not match the VCF which was actually written out
    """
    # Read in the VCF (that was just written out) back in with pysam
    flagged_contents = VariantFile(flagged_vcf, 'r')

    # Get list of lines in original header that was written out
    # Get list of variant records which were to be written out
    # Make one big list for original VCF contents which were to be written out
    original_header = list(set(
        str(header) for header in original_vcf_contents.header.records
    ))
    original_records = [
        str(var) for variant_list in gene_variant_dict.values()
        for var in variant_list
    ]
    original_contents = original_header + original_records

    # Get list of lines which were written out to VCF file header
    # Get list of variant records which were in the written out VCF
    # Make one big list for written out VCF
    written_header = list(set(
        str(header) for header in flagged_contents.header.records
    ))
    written_records = [
        str(variant_record) for variant_record in flagged_contents
    ]
    written_contents = written_header + written_records

    assert original_contents == written_contents, (
        "Header and variants written to VCF not identical to those "
        "intended to be written out"
    )


def bcftools_sort(input_vcf):
    """
    Sort the VCF. As we write out variants in gene order, if a
    variant is annotated to multiple genes/transcripts and split, we can end up
    with the VCF not being sorted by position so this fixes that

    Parameters
    ----------
    input_vcf : str
        path to filtered VCF to sort
    output_vcf : str
        name of output sorted VCF

    Outputs
    -------
    {vcf}.sorted.vcf : file
        sorted vcf file output from bcftools
    """
    print(f"Sorting flagged VCF {input_vcf} with bcftools sort")
    output_vcf = f"{Path(input_vcf).stem.split('.')[0]}.sorted.vcf"

    # Check total rows before running bcftools sort
    pre_sort_total = subprocess.run(
        f"zgrep -v '^#' {input_vcf} | wc -l",
        shell=True,
        capture_output=True
    )

    cmd = f"bcftools sort {input_vcf} -o {output_vcf}"

    output = subprocess.run(cmd, shell=True, capture_output=True)

    assert output.returncode == 0, (
        f"\n\tError in sorting VCF with bcftools. VCF: {input_vcf}"
        f"\n\tExitcode:{output.returncode}"
        f"\n\t{output.stderr.decode()}"
    )

    # Check total rows after sorting
    post_sort_total = subprocess.run(
        f"zgrep -v '^#' {output_vcf} | wc -l",
        shell=True,
        capture_output=True
    )

    pre_sort = pre_sort_total.stdout.decode()
    post_sort = post_sort_total.stdout.decode()

    print(
        f"Total lines before sorting: {pre_sort}"
        f"Total lines after sorting: {post_sort}"
    )

    assert pre_sort == post_sort, (
        "Count of variants before and after bcftools sort do not match"
    )

def bgzip(file) -> None:
    """
    Call bgzip on a given file

    Parameters
    ----------
    file : file to compress

    Outputs
    -------
    compressed file

    Raises
    ------
    AssertionError
        Raised when non-zero exit code returned by bgzip
    """
    print(f"Calling bgzip on {file}")

    output = subprocess.run(
        f"bgzip --force {file} -c > {file}.gz",
        shell=True,
        capture_output=True
    )

    assert output.returncode == 0, (
        f"\n\tError in compressing file with bgzip. File: {file}"
        f"\n\tExitcode:{output.returncode}"
        f"\n\t{output.stderr.decode()}"
    )



def add_annotation(input_file, transc):
    """
    Main function to take a VCF and add the INFO field required for filtering

    Parameters
    ----------
    input_vcf : str
        name of the input VCF
    panel_dict : dict
        default dict with each gene on panel as key and the gene info as val
    filter_command : str
        full bcftools filter command
    """
    annotated_vcf = f"{Path(input_file)}"
    split_vcf = f"{Path(input_file).stem.split('.')[0]}.split.vcf"
    MANE_flagged_vcf = f"{Path(input_file).stem.split('.')[0]}.flagged.vcf"
    sorted_vcf = f"{Path(input_file).stem.split('.')[0]}.sorted.vcf"

    #decompress .vcf.gz file from vep output
    input_vcf = decompress(input_file)
  
    # separate csq fields (creates split_vcf)
    bcftools_pre_process(input_vcf)

    # create pysam object of vcf for flagging
    vcf_contents, sample_name = read_in_vcf(split_vcf)

    # add MANE field from config
    gene_var_dict = add_MANE_field(vcf_contents, transc)
    write_out_flagged_vcf(MANE_flagged_vcf, gene_var_dict, vcf_contents)
    check_written_out_vcf(vcf_contents, gene_var_dict, MANE_flagged_vcf)
    bcftools_sort(MANE_flagged_vcf)

    bgzip(MANE_flagged_vcf)
    os.remove(annotated_vcf)
    os.remove(split_vcf)
    os.remove(sorted_vcf)




