#!/bin/bash
#SBATCH --job-name=counts
#SBATCH --cpus-per-task=50
#SBATCH --mem=200G
#SBATCH --time=23:00:00
#SBATCH --output=$(pwd)/counts_%j.out
#SBATCH --error=$(pwd)/counts_%j.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=USERNAME@student.le.ac.uk
#SBATCH --export=NONE

# Get the current working directory
work_dir=$(pwd)

module load R/4.3.1

R --vanilla <<EOF

# Check and install required packages if not already installed
if (!requireNamespace("ensembldb", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("ensembldb")
}

if (!requireNamespace("Rsubread", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("Rsubread")
}

# Load required packages
library(Rsubread)
library(ensembldb)

# Set input and output paths relative to the current working directory
bam_dir <- file.path("$work_dir", "bam_files")
output_dir <- file.path("$work_dir", "featureCounts_output")
ref_dir <- "$work_dir"

# Get sorted BAM files
bam_files <- list.files(path=bam_dir, pattern="*_sorted.bam$", full.names=TRUE)

# Download and prepare the GTF file
gtf_file <- file.path(ref_dir, "Homo_sapiens.GRCh38.96.gtf")
# Download the GTF file and place it in the ref_dir

# Run featureCounts
bamcounts <- featureCounts(
  files = bam_files,
  annot.ext = gtf_file,
  isGTFAnnotationFile = TRUE,
  GTF.featureType = "exon",
  GTF.attrType = "gene_id",
  useMetaFeatures = TRUE,
  allowMultiOverlap = TRUE,
  isPairedEnd = TRUE,
  requireBothEndsMapped = TRUE,
  countChimericFragments = FALSE,
  nthreads = 30
)

# Write featureCounts output to files in the output directory
for (n in names(bamcounts)){
  write.table(bamcounts[[n]], file=file.path(output_dir, paste0(n, ".csv")), sep=",", quote=F, col.names=NA)
}
EOF
