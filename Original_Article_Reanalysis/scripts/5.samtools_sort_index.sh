#!/bin/bash

# Load the Samtools module
module load samtools/1.17-wenuvv5

# Define directories
ALIGN_DIR="$(pwd)/STAR_alignment/aligned_outputs"
BAM_DIR="$(pwd)/STAR_alignment/bam_files"

# Create the output directory
mkdir -p $BAM_DIR

# Set the number of threads for Samtools (locally there are 64 cpus)
THREADS=60

# Process ALL SAM files
for sam in $ALIGN_DIR/*_Aligned.out.sam; do
    base=$(basename $sam _Aligned.out.sam)
    samtools view -bS $sam | samtools sort -@ $THREADS -o $BAM_DIR/${base}_sorted.bam
    samtools index $BAM_DIR/${base}_sorted.bam
done
