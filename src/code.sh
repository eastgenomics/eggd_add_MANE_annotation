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
    python3 /add_MANE_annotation.py \
        -i $input_file_path \
        -f $transcript_file_path 
    # prepare outputs
    echo "All scripts finished successfully, uploading output files to dx"
    if [ "$debug_fail_end" == 'true' ]; then exit 1; fi
    mkdir -p ~/out/result_files/
    mv *.vcf.gz ~/out/result_files/

    # Upload output files
    dx-upload-all-outputs --parallel

}
