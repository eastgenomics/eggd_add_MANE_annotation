#!/bin/bash
# Add_MANE_annotation
#
# Any code outside of main() (or any entry point you may add) is
# ALWAYS executed, followed by running the entry point itself.

# Exit at any point if there is any error and output each line as it is executed (for debugging)
set -e -x -o pipefail

# install dependencies
sudo -H python3 -m pip install --no-index --no-deps packages/*
export BCFTOOLS_PLUGINS=/usr/local/libexec/bcftools


main() {

    # get inputs
    dx-download-all-inputs --parallel

    # run tool
    python3 add_MANE_annotation.py \
        -i $input_vcf_annotated_path \
        -t $transcript_file_path 
    # prepare outputs
    echo "All scripts finished successfully, uploading output files to dx"
    mkdir -p ~/out/output_vcf_annotated_mane/
    mv *.vcf.gz ~/out/output_vcf_annotated_mane/

    # Upload output files
    dx-upload-all-outputs --parallel

}
