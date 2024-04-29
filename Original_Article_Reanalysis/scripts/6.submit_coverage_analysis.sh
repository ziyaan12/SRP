#!/bin/bash

# Define directories
BAM_DIR="$(pwd)/bam_files"
COVERAGE_DIR="$(pwd)/coverage_tracks"
COUNTS_DIR="$(pwd)/read_counts"

# Get list of sorted BAM files
BAM_FILES=($BAM_DIR/*sorted.bam)

# Define the number of jobs to create
TOTAL_FILES=${#BAM_FILES[@]}
FILES_PER_JOB=$((TOTAL_FILES / 80))
COUNTER=0
JOB_NUMBER=1

# Create separate folders for job outputs to keep things organized
mkdir -p $(pwd)/coverage_analysis_jobs

# Generate job scripts for each subset of BAM files
for ((i=0; i<TOTAL_FILES; i+=FILES_PER_JOB)); do
    # Define a batch script for each subset
    BATCH_SCRIPT=$(pwd)/coverage_analysis_jobs/coverage_job_$JOB_NUMBER.sbatch
    
    # Generate the script
    echo "#!/bin/bash" > $BATCH_SCRIPT
    echo "#SBATCH --job-name=coverage_analysis_$JOB_NUMBER" >> $BATCH_SCRIPT
    echo "#SBATCH --cpus-per-task=10" >> $BATCH_SCRIPT
    echo "#SBATCH --mem=200G" >> $BATCH_SCRIPT
    echo "#SBATCH --time=5:00:00" >> $BATCH_SCRIPT
    echo "#SBATCH --output=$(pwd)/coverage_analysis_jobs/coverage_$JOB_NUMBER_%j.out" >> $BATCH_SCRIPT
    echo "#SBATCH --error=$(pwd)/coverage_analysis_jobs/coverage_$JOB_NUMBER_%j.err" >> $BATCH_SCRIPT
    echo "#SBATCH –mail-type=END,FAIL" >> $BATCH_SCRIPT
    echo "#SBATCH --mail-user=USERNAME@student.le.ac.uk" >> $BATCH_SCRIPT
    echo "#SBATCH --export=NONE" >> $BATCH_SCRIPT
    echo "" >> $BATCH_SCRIPT
    echo "module load R/4.3.1" >> $BATCH_SCRIPT
    echo "" >> $BATCH_SCRIPT
    
    echo "R --vanilla <<EOF" >> $BATCH_SCRIPT
    echo "# Load required packages" >> $BATCH_SCRIPT
    echo "library(bamsignals)" >> $BATCH_SCRIPT
    echo "library(rtracklayer)" >> $BATCH_SCRIPT
    echo "library(Rsubread)" >> $BATCH_SCRIPT
    echo "library(GenomicRanges)" >> $BATCH_SCRIPT
    echo "library(Rsamtools)" >> $BATCH_SCRIPT
    echo "" >> $BATCH_SCRIPT
    
    for ((j=i; j<i+FILES_PER_JOB && j<TOTAL_FILES; j++)); do
        bam=${BAM_FILES[j]}
        sample=$(basename $bam _sorted.bam)
        
        echo "# Process sample: $sample" >> $BATCH_SCRIPT
        echo "bam <- \"$bam\"" >> $BATCH_SCRIPT
        echo "sample <- sub(\"_sorted.bam\", \"\", basename(bam))" >> $BATCH_SCRIPT
        echo "" >> $BATCH_SCRIPT
        echo "# Extract sequence information from BAM header" >> $BATCH_SCRIPT
        echo "seq_info <- seqinfo(BamFile(bam))" >> $BATCH_SCRIPT
        echo "" >> $BATCH_SCRIPT
        echo "# Define the genome tiles" >> $BATCH_SCRIPT
        echo "tiles <- tileGenome(seq_info, tilewidth = 20, cut.last.tile.in.chrom = TRUE)" >> $BATCH_SCRIPT
        echo "" >> $BATCH_SCRIPT
        echo "# Count 5' ends in 20bp windows" >> $BATCH_SCRIPT
        echo "cov <- bamCount(bam, tiles, shift = 0, ss = TRUE)" >> $BATCH_SCRIPT
        echo "" >> $BATCH_SCRIPT
        echo "# Convert coverage counts to GRanges objects" >> $BATCH_SCRIPT
        echo "cov_fwd <- tiles" >> $BATCH_SCRIPT
        echo "score(cov_fwd) <- cov[\"sense\", ]" >> $BATCH_SCRIPT
        echo "strand(cov_fwd) <- \"+\"" >> $BATCH_SCRIPT
        echo "" >> $BATCH_SCRIPT
        echo "cov_rev <- tiles" >> $BATCH_SCRIPT
        echo "score(cov_rev) <- cov[\"antisense\", ]" >> $BATCH_SCRIPT
        echo "strand(cov_rev) <- \"-\"" >> $BATCH_SCRIPT
        echo "" >> $BATCH_SCRIPT
        echo "# Write bigWig files" >> $BATCH_SCRIPT
        echo "export.bw(cov_fwd, file.path(\"$COVERAGE_DIR\", paste0(sample, \"_fwd.bw\")))" >> $BATCH_SCRIPT
        echo "export.bw(cov_rev, file.path(\"$COVERAGE_DIR\", paste0(sample, \"_rev.bw\")))" >> $BATCH_SCRIPT
        echo "" >> $BATCH_SCRIPT
    done
    
    echo "EOF" >> $BATCH_SCRIPT
    
    chmod +x $BATCH_SCRIPT
    
    # Increment the job number
    ((JOB_NUMBER++))
done

# Submit all jobs in order
for script in $(ls -v $(pwd)/coverage_analysis_jobs/*.sbatch); do
    sbatch $script
done
