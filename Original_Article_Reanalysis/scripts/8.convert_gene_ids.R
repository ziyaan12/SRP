# Load necessary libraries for annotation
library(AnnotationDbi)
library(EnsDb.Hsapiens.v79)

# Read count data from CSV and extract Ensembl IDs
count_data <- read.csv(file = 'counts.csv', row.names = 1)
ensembl_ids <- rownames(count_data)

# Use ensembldb to fetch corresponding gene symbols for the Ensembl IDs
gene_info <- ensembldb::select(EnsDb.Hsapiens.v79, keys=ensembl_ids, 
                               keytype='GENEID', columns='SYMBOL')

# Construct a named vector for gene symbols using Ensembl IDs as the names
symbol_vector <- setNames(gene_info$SYMBOL, gene_info$GENEID)

# Align gene symbols with Ensembl IDs in the count data
aligned_symbols <- symbol_vector[rownames(count_data)]

# Append Ensembl IDs to non-unique gene symbols to ensure uniqueness
non_unique <- duplicated(aligned_symbols) | duplicated(aligned_symbols, fromLast = TRUE)
aligned_symbols[non_unique] <- paste(aligned_symbols[non_unique], rownames(count_data)[non_unique], sep="_")

# Update row names in the count data with the aligned, unique gene symbols
rownames(count_data) <- aligned_symbols

# Save the modified count data to a new CSV file
write.csv(count_data, "counts_updated.csv", row.names = TRUE)
