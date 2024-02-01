eggd_add_MANE_annotation
This app takes an annotated VCF (output of eggd_VEP) and adds MANE info to VCF

Usage
To run the app:

dx run app-GZ9FZ78457v7qjBXPXqGByyP \
    -iinput_vcf=[annotated vcf] \
    -itranscript_file=[MANE file] \
\
    --destination=/path/to/output/dir -y

