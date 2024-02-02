"""
Main script which takes a file and adds our flag for optimised filtering
"""
import argparse
import re
import pandas as pd

import vcf


def parse_args() -> argparse.Namespace:
    """
    Parse the command line arguments inputs given

    Returns
    -------
    args : Namespace
        Namespace of passed command line argument inputs
    """

    parser = argparse.ArgumentParser(
        description='Information necessary for add MANE annotation'
    )

    # Add CLI arg of the input annotated VCF to have filtering flags added
    parser.add_argument(
        '-i',
        '--input_vcf_annotated',
        type=str,
        required=True,
        help='Annotated VCF file to have MANE annotation added'
    )

    # Add CLI input of MANE transcripts
    parser.add_argument(
        '-t',
        '--transcript_file',
        type=str,
        required=True,
        help='File with MANE transcripts to be converted into a list'
    )

    args = parser.parse_args()

    return args

def get_list_transcripts(transcript_file):
    """
    read the MANE file table and convert the transcripts column into a list

    parameters
    ----------
    transcript_file: file of MANE transcripts in tsv format

    output
    ----------
    list of mane transcripts    
    """
    table_mane = pd.read_csv(transcript_file, sep="\t")
    transcript_list = table_mane['RefSeq_nuc'].to_list()

    return transcript_list

def main():
    args = parse_args()
    input_vcf_decompressed = vcf.decompress(
        args.input_vcf_annotated
    )
    transcript_list = get_list_transcripts(args.transcript_file)
    vcf.add_annotation(
        input_vcf_decompressed, transcript_list
    )
    vcf.compress(
        args.MANE_flagged_vcf
    )

if __name__ == "__main__":
    main()