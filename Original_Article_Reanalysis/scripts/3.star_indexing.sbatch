#!/bin/bash
#SBATCH --job-name=star_indexing
#SBATCH --cpus-per-task=30
#SBATCH --mem=40G
#SBATCH --time=01:00:00
#SBATCH --output=star_index_%j.out
#SBATCH --error=star_index_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=USERNAME@student.le.ac.uk
#SBATCH --export=NONE

#Load the STAR module
module load star/2.7.10b-m3zkpic

# Define directories
REFERENCE_DIR="$(pwd)/STAR_alignment/reference"
INDEX_DIR="$(pwd)/STAR_alignment/STAR_indices"
ALIGN_DIR="$(pwd)/STAR_alignment/aligned_outputs"
FASTQ_DIR="$(pwd)/SRP"

# Generate STAR genome indices
STAR --runMode genomeGenerate \
     --runThreadN 30 \
     --genomeDir $INDEX_DIR \
     --genomeFastaFiles $REFERENCE_DIR/hg38.fa \
     --genomeSAindexNbases 14
