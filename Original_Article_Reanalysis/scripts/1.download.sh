#!/bin/bash

# Loop through each accession number in the file
while read -r SRR; do
    echo "Processing: $SRR"
    
    # Download the .sra file
    echo "Downloading .sra file for: $SRR"
    prefetch $SRR

    # Convert the .sra file to FASTQ format
    echo "Converting to FASTQ: $SRR"
    fasterq-dump $SRR
    
    # Remove the directory and .sra file
    echo "Cleaning up: $SRR"
    rm -rf $SRR
    
done < SRR_Acc_List.txt

echo "All downloads and conversions are complete."
