# eggd_add_MANE_annotation
This app takes an annotated VCF (output of eggd_VEP) and adds MANE info to VCF 
## Usage

To run the app:

```
dx run app-GZ9FZ78457v7qjBXPXqGByyP \
    -i input_vcf_annotated \
    -t transcript_file \

```

# Local Tool README

## add_MANE_annotation
Add MANE annotation is a tool used to add MANE info.

### Description
INFO field will be added (named 'MANE')


Add MANE annotation uses:
- [bcftools](https://samtools.github.io/bcftools/bcftools.html, "bcftools website")
- [pysam](https://pysam.readthedocs.io/en/latest/, "pysam documentation")

### Requirements
- bcftools
    - bcftools +split-vep

