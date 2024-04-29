cd $SCRATCHDIR/STAR_alignment/fastq-screen

echo "DATABASE hg38 /scratch/alice/U/USERNAME/STAR_alignment/fastq-screen/hg38
DATABASE mm10 /scratch/alice/U/USERNAME/STAR_alignment/fastq-screen/mm10" > fastq_screen.conf

echo "module load bowtie2/2.5.1-s4zfmrn

for f in /scratch/alice/b/bmsas1/SRP/*_1.fastq; do
    fastq_screen --aligner bowtie2 --conf /scratch/alice/U/USERNAME/STAR_alignment/fastq-screen/fastq_screen.conf --outdir /scratch/alice/U/USERNAME/STAR_alignment/fastq-screen/fastq_screen_out \$f \${f/_1.fastq/_2.fastq}
done" > fastq_screen.bash
