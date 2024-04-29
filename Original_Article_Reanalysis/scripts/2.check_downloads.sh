#!/bin/bash

# Loop through each accession number in the file
while read -r SRR; do
    # Check if both paired-end FASTQ files exist
    if [ ! -f "${SRR}_1.fastq" ] || [ ! -f "${SRR}_2.fastq" ]; then
        # If either file is missing, print the SRR accession number
        echo "Missing files for: $SRR"
    fi
done < SRR_Acc_List.txt
echo "Check complete."
