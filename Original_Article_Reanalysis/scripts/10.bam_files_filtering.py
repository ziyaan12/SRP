import os
import shutil

def get_human_mouse_ratio(fastq_screen_file):
    with open(fastq_screen_file, 'r') as f:
        human_mapped = 0
        mouse_mapped = 0
        for line in f:
            if line.startswith('hg38'):
                parts = line.split()
                human_mapped = float(parts[5]) + float(parts[7])  # Sum of %One_hit_one_genome and %Multiple_hits_one_genome for human
            elif line.startswith('mm10'):
                parts = line.split()
                mouse_mapped = float(parts[5]) + float(parts[7])  # Sum of %One_hit_one_genome and %Multiple_hits_one_genome for mouse
    return human_mapped / mouse_mapped if mouse_mapped > 0 else float('inf')

bam_dir = '/scratch/alice/U/USERNAME/STAR_alignment/bam_files/'
filtered_bam_dir = '/scratch/alice/U/USERNAME/STAR_alignment/filtered_bam_files/'

# Create directory for filtered bam files
if not os.path.exists(filtered_bam_dir):
    os.makedirs(filtered_bam_dir)


# Copy all bam and bai files to the new directory
for file_name in os.listdir(bam_dir):
    if file_name.endswith('.bam') or file_name.endswith('.bai'):
        shutil.copy(os.path.join(bam_dir, file_name), os.path.join(filtered_bam_dir, file_name))

# Find bam files to remove based on the human-to-mouse ratio
bam_files_to_remove = set()
for fastq_screen_file in os.listdir('fastq_screen_out'):
    if fastq_screen_file.endswith('.txt'):
        ratio = get_human_mouse_ratio(f'fastq_screen_out/{fastq_screen_file}')
        if ratio < 5:
            sample = fastq_screen_file.split('_')[0]
            bam_file = os.path.join(filtered_bam_dir, f'{sample}_sorted.bam')
            if os.path.exists(bam_file):
                bam_files_to_remove.add(bam_file)

# Remove specified bam and corresponding bai files
for bam_file in bam_files_to_remove:
    os.remove(bam_file)
    bai_file = bam_file + '.bai'
    if os.path.exists(bai_file):
        os.remove(bai_file)

# Print results
bam_files_to_remove = sorted(bam_files_to_remove)  # Convert set to sorted list
print('\n'.join(bam_files_to_remove))
print(f'filtered_bam_files directory was created successfully. Total number of BAM files removed: {len(bam_files_to_remove)}')
