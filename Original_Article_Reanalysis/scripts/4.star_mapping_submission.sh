#!/bin/bash

# Define directories
REFERENCE_DIR="$(pwd)/STAR_alignment/reference"
INDEX_DIR="$(pwd)/STAR_alignment/STAR_indices"
ALIGN_DIR="$(pwd)/STAR_alignment/aligned_outputs"
FASTQ_DIR="$(pwd)/SRP"
FASTQ_FILES=($FASTQ_DIR/*_1.fastq)

# Define how many jobs you want to create (total number of FASTQ pairs divided by the number you want per job).
TOTAL_FILES=${#FASTQ_FILES[@]}
FILES_PER_JOB=$((TOTAL_FILES / 15))
COUNTER=0
JOB_NUMBER=1

# Create separate folders for job outputs to keep things organized.
mkdir -p $(pwd)/star_mapping_jobs

for ((i=0; i<TOTAL_FILES; i+=FILES_PER_JOB)); do
    # Define a batch script for each subset.
    BATCH_SCRIPT=$(pwd)/star_mapping_jobs/star_job_$JOB_NUMBER.sbatch
    
    # Generate the script.
    echo "#!/bin/bash" > $BATCH_SCRIPT
    echo "#SBATCH --job-name=star_alignment_$JOB_NUMBER" >> $BATCH_SCRIPT
    echo "#SBATCH --cpus-per-task=30" >> $BATCH_SCRIPT
    echo "#SBATCH --mem=60G" >> $BATCH_SCRIPT
    echo "#SBATCH --time=13:00:00" >> $BATCH_SCRIPT
    echo "#SBATCH --output=$(pwd)/star_mapping_jobs/star_$JOB_NUMBER_%j.out" >> $BATCH_SCRIPT
    echo "#SBATCH --error=$(pwd)/star_mapping_jobs/star_$JOB_NUMBER_%j.err" >> $BATCH_SCRIPT
    echo "#SBATCH --mail-type=END,FAIL" >> $BATCH_SCRIPT
    echo "#SBATCH --mail-user=USERNAME@student.le.ac.uk" >> $BATCH_SCRIPT
    echo "#SBATCH --export=NONE" >> $BATCH_SCRIPT
    echo "" >> $BATCH_SCRIPT
    echo "module load star/2.7.10b-m3zkpic" >> $BATCH_SCRIPT
    echo "" >> $BATCH_SCRIPT

    # Add commands for STAR alignment for each subset of files.
    for ((j=i; j<i+FILES_PER_JOB; j++)); do
        fastq1=${FASTQ_FILES[j]}
        fastq2="${fastq1%_1.fastq}_2.fastq"
        sample_id=$(basename $fastq1 _1.fastq)
        
        echo "STAR --runMode alignReads \\" >> $BATCH_SCRIPT
        echo "     --runThreadN 30 \\" >> $BATCH_SCRIPT
        echo "     --genomeDir $INDEX_DIR \\" >> $BATCH_SCRIPT
        echo "     --readFilesIn $fastq1 $fastq2 \\" >> $BATCH_SCRIPT
        echo "     --sjdbGTFfile $REFERENCE_DIR/genes.gtf \\" >> $BATCH_SCRIPT
        echo "     --sjdbOverhang 100 \\" >> $BATCH_SCRIPT
        echo "     --outFileNamePrefix $ALIGN_DIR/${sample_id}_ \\" >> $BATCH_SCRIPT
        echo "     --outSAMstrandField intronMotif \\" >> $BATCH_SCRIPT
        echo "     --outFilterIntronMotifs RemoveNoncanonical" >> $BATCH_SCRIPT
        echo "" >> $BATCH_SCRIPT
    done

    chmod +x $BATCH_SCRIPT

    # Increment the job number.
    ((JOB_NUMBER++))
done

# Submit all jobs
for script in $(pwd)/star_mapping_jobs/*.sbatch; do
    sbatch $script
done
