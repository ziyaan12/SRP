
![Logo](https://img.hotimg.com/Screenshot-2024-04-24-121024.md.jpeg)


# Project Title
Re-analysis of Metastatic Breast Cancer scRNA-seq Data




## API Reference
https://github.com/ziyaan12/SRP




## Appendix

    Project Overview:

This project aims to re-analyze the single-cell RNA sequencing (scRNA-seq) data from the study "Single-cell Analysis Reveals Inter- and Intratumour Heterogeneity in Metastatic Breast Cancer" by Hamelin et al. (2023). The primary objectives are:

1. To gain a deeper understanding of the molecular and cellular heterogeneity within metastatic breast cancer tumors.
2. To identify novel therapeutic targets or prognostic biomarkers that could improve the clinical management of metastatic breast cancer.
3. To develop a user-friendly web interface that allows researchers, clinicians, and the broader scientific community to explore and interact with the re-analyzed dataset.

Understanding the complex landscape of metastatic breast cancer at the single-cell level is crucial for developing more effective personalized treatment strategies. By leveraging the rich dataset from the Hamelin et al. study, this project will provide new insights into the molecular drivers of metastasis and tumor evolution, potentially leading to improved patient outcomes.

    Data Sources:

1. **Original Dataset**: The scRNA-seq data from the Hamelin et al. study is available in the Gene Expression Omnibus (GEO) database under the accession number [GSE202695](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE202695).
2. **Raw Sequencing Data**: The raw RNA sequencing data can be accessed from the NCBI Short Read Archive (SRA) under the BioProject accession [PRJNA837007](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA837007).
3. **Reference Genome and Annotation**: The human reference genome (hg38) and gene annotation file (genes.gtf) will be used for the data processing and analysis steps.



       Software Requirements:

1.GEO (Gene Expression Omnibus) website: Accessible through a web browser such as Firefox.

2.SRA (Sequence Read Archive) Run Selector: Accessible through a web browser such as Firefox.

3.NoMachine: Remote desktop software used to connect to ALICE HPC.

4.Terminal: Command-line interface for executing commands on ALICE HPC.

5.SRA Toolkit: A collection of tools for working with SRA data. It should be pre-installed on ALICE HPC.

6.STAR (Spliced Transcripts Alignment to a Reference) software: A computational tool for performing read alignment. It should be pre-installed on ALICE HPC.

7.Text editor: Used for creating and editing bash scripts.


    Usage Instructions:

1. **Download Data from GEO (GSE202695)**: Navigate to the GEO homepage, search for the dataset with accession number GSE202695, and download the available processed data files.

2. **Download Raw Data from SRA (PRJNA837007) on ALICE HPC**: Log into the ALICE HPC system, open Firefox, and navigate to the SRA Run Selector to download the accession list file (SRR_Acc_List.txt). Then, create a bash script to download and convert the SRA files to FASTQ format.

3. **STAR Indexing on ALICE HPC**: Create a batch script to generate the STAR genome indices using the reference genome (hg38.fa) and annotation file (genes.gtf).

4. **STAR Mapping on ALICE HPC**: Create a bash script that generates multiple job submission scripts to perform the STAR alignment in parallel.

5. **Sorting and Indexing using Samtools on ALICE HPC**: Create a bash script to sort and index the SAM files using Samtools.

6. **Generating Strand-Specific Coverage Tracks Using R on ALICE HPC**: Create an R script to generate strand-specific coverage tracks from the BAM files and write them as bigWig files. Submit the R script as a batch job on ALICE HPC.


    File Structure:   


   
1.**Scratch Directory**: This is the main working directory, which can be accessed using the `$SCRATCHDIR` environment variable.

2.Inside the Scratch Directory:
   - `SRP`: This directory is created to store the downloaded raw data files (SRA/FASTQ) from the SRA.
   - `STAR_alignment`: This is the main directory for organizing the data processing steps. It contains the following subdirectories:
     - `reference`: This directory stores the reference genome (hg38.fa) and the gene annotation file (genes.gtf).
     - `STAR_indices`: This directory will store the STAR genome indices generated during the indexing step.
     - `aligned_outputs`: This directory will store the STAR-aligned SAM/BAM files.
     - `bam_files`: This directory will store the sorted and indexed BAM files after processing with Samtools.
     - `coverage_tracks`: This directory will store the strand-specific coverage tracks generated using the R script.
     - `read_counts`: This directory will store the read count data generated during the coverage track creation.
     - `star_mapping_jobs`: This directory will contain the individual job submission scripts for parallel STAR alignment.



 

## Authors



Bashar M.S. Al-Smadi,
Amirhossein Soltani,
Shyam S. Gopi Krishnan,
Esther Jacob,
Dedeepya Reddy Marthala,
Noha Nashnoush,
Ziyaan Osman,
Sanaa Patel

## License

[MIT](https://choosealicense.com/licenses/mit/)
[APACHE](https://www.apache.org/licenses/LICENSE-2.0)


## Badges



[![MIT License](https://img.shields.io/badge/License-MIT-green.svg)](https://choosealicense.com/licenses/mit/)
[![GPLv3 License](https://img.shields.io/badge/License-GPL%20v3-yellow.svg)](https://opensource.org/licenses/)
[![AGPL License](https://img.shields.io/badge/license-AGPL-blue.svg)](http://www.gnu.org/licenses/agpl-3.0)


## References

Hamelin, B., Milan, Sethi, A., Kloc, M., Münst, S., Beisel, C., Eschbach, K., Kohler, H., Savas Soysal, Vetter, M., Weber, W.P., Stadler, M.B. and Bentires-Alj, M. (2023). Single-cell Analysis Reveals Inter- and Intratumour Heterogeneity in Metastatic Breast Cancer. Journal of mammary gland biology and neoplasia, 28(1). doi:https://doi.org/10.1007/s10911-023-09551-z.

‌