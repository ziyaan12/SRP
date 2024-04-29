# Install devtools if not already installed
if (!require("devtools")) install.packages("devtools")
# Install the QUIC package from source
devtools::install_url("https://cran.r-project.org/src/contrib/Archive/QUIC/QUIC_1.1.tar.gz", repos = NULL, type = "source")
# Install the griph package from GitHub
devtools::install_github("pappasakis/griph/griph")

# Load required libraries
library(data.table)
library(griph)
library(scater)
library(AnnotationDbi)
library(stats)
library(org.Hs.eg.db)

# Set working directory (Session > Set Working Directory > To Source File Location)

# Load the filtered counts data
counts <- fread("re_analysis_counts_afterQC.csv")

# Remove rows with NaN and duplicates
counts <- na.omit(counts)  
counts <- unique(counts)

# Convert gene symbols to Entrez IDs
human <- org.Hs.eg.db
entrez_ids <- mapIds(human, keys = counts$external_gene_name, keytype = "SYMBOL", column = "ENTREZID")
counts$external_gene_name <- entrez_ids

# Extract the count matrix and set row names to Entrez IDs
count_matrix <- as.matrix(counts[, -1]) 
rownames(count_matrix) <- counts$external_gene_name

# Log-normalize the counts
log_counts <- log2(count_matrix + 1)

# Infer cell cycle stages
cc_scores <- predictCellCycle(log_counts, org = "human.Whitfield", cor_thr = 0.2, refine_iter = 200)

# Calculate library complexity (percentage of detected genes)
lib_complexity <- colSums(log_counts > 0) / nrow(log_counts) * 100

# Create a design matrix with cell cycle scores and library complexity
design_matrix <- model.matrix(~ cc_scores + lib_complexity)

# Regress out cell cycle and library complexity using lm
reg_counts <- t(apply(log_counts, 1, function(y) {
  fit <- lm(y ~ design_matrix)
  residuals(fit)
}))

# Add row and column names to the corrected counts matrix
rownames(reg_counts) <- rownames(log_counts)
colnames(reg_counts) <- colnames(log_counts)

# Convert Entrez IDs back to Ensembl gene names
ensembl_ids <- mapIds(human, keys = rownames(reg_counts), keytype = "ENTREZID", column = "SYMBOL")
rownames(reg_counts) <- ensembl_ids

# Save the cell cycle corrected log counts to a CSV file
fwrite(as.data.table(reg_counts, keep.rownames = "external_gene_name"), 
       "re_analysis_CCcorrected_logcounts.csv")
