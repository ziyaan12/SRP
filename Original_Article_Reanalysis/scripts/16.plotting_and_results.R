# Install devtools if not already installed
if (!require("devtools")) install.packages("devtools")

# Install the QUIC package from source
devtools::install_url("https://cran.r-project.org/src/contrib/Archive/QUIC/QUIC_1.1.tar.gz", repos = NULL, type = "source")

# Install the griph package from GitHub
devtools::install_github("ppapasaikas/griph/griph")

library(griph)

install.packages("BiocManager")
BiocManager::install("scater")
BiocManager::install("scran")
BiocManager::install("fgsea")
BiocManager::install("msigdb")
BiocManager::install("msigdbr")

# Load required packages
library(scater)
library(scran)
library(Rtsne)
library(umap)
library(igraph)
library(plotly)



# Read in the log-normalized counts data
counts <- read.csv("re_analysis_CCcorrected_logcounts.csv", row.names = 1)

# Create a SingleCellExperiment object
sce <- SingleCellExperiment(assays = list(logcounts = as.matrix(counts)))


gc()

# Compute PCA
sce <- runPCA(sce, exprs_values = "logcounts")

# Compute t-SNE
set.seed(123)
sce <- runTSNE(sce, dimred = "PCA", perplexity = 30)

# Compute UMAP
set.seed(123)
sce <- runUMAP(sce, dimred = "PCA")

# Build nearest-neighbor graph using t-SNE
g_tsne <- buildSNNGraph(sce, use.dimred = "TSNE")

# Extract model IDs from cell names
#header row format: `s<ID>_mouse<ID>_<modelID>_<tissueType>_<cellDescription>`
sce$model <- ifelse(grepl("^s\\d+_mouse\\d+_HBRX\\d+_", colnames(sce)),
                    sub("^s\\d+_mouse\\d+_(HBRX\\d+)_.*", "\\1", colnames(sce)),
                    ifelse(grepl("^s\\d+_mouse\\d+_MDAMB231_", colnames(sce)),
                           "MDAMB231", "MDAMB231"))

sce$model <- sub("HBRX3078", "PDX1", sce$model)
sce$model <- sub("HBRX2353", "PDX2", sce$model)
sce$model <- sub("HBRX1921", "PDX3", sce$model)
sce$model <- sub("HBRX2344", "PDX4", sce$model)

# Extract tissue types from cell names
sce$tissue <- sub("^s\\d+_mouse\\d+_\\w+_(.*?)_\\w+cells$", "\\1", colnames(sce))

# Normalize tissue types
sce$tissue <- gsub("Lung.CTRL", "Lung", sce$tissue)
sce$tissue <- gsub("Tumor.Control", "Tumor", sce$tissue)
sce$tissue <- gsub("Metastases", "Lung", sce$tissue)
sce$tissue <- gsub("s0024_mouse1a_MDAMB231_Tumor_1cells", "Tumor", sce$tissue)
sce$tissue <- gsub("s0035_mouse1a_MDAMB231_Tumor_1cells", "Tumor", sce$tissue)
sce$tissue <- gsub("_\\w+$", "", sce$tissue)


# Plot t-SNE colored by model (Fig. 1d equivalent)
plotTSNE(sce, colour_by = "model") + ggtitle("t-SNE by Model")

# Plot t-SNE colored by tissue (Fig. 1e equivalent) 
plotTSNE(sce, colour_by = "tissue") + ggtitle("t-SNE by Tissue")

# Cluster cells using Louvain algorithm on t-SNE graph
set.seed(123)
sce$clusters_tsne <- factor(cluster_louvain(g_tsne)$membership)

# Plot t-SNE colored by clusters (Fig. 2c equivalent)
plotTSNE(sce, colour_by = "clusters_tsne") + ggtitle("t-SNE Clusters")


#¬
# Plot interactive t-SNE colored by model (Fig. 1d equivalent)
tsne_model <- plotTSNE(sce, colour_by = "model")
tsne_model_interactive <- ggplotly(tsne_model, tooltip = "text") %>%
  layout(title = "t-SNE by Model")

# Print the interactive plot for model
tsne_model_interactive

# Plot interactive t-SNE colored by tissue (Fig. 1e equivalent)
tsne_tissue <- plotTSNE(sce, colour_by = "tissue")
tsne_tissue_interactive <- ggplotly(tsne_tissue, tooltip = "text") %>%
  layout(title = "t-SNE by Tissue")

# Print the interactive plot for tissue
tsne_tissue_interactive

# Plot interactive t-SNE colored by clusters (Fig. 2c equivalent)
tsne_clusters <- plotTSNE(sce, colour_by = "clusters_tsne")
tsne_clusters_interactive <- ggplotly(tsne_clusters, tooltip = "text") %>%
  layout(title = "t-SNE Clusters")

# Print the interactive plot for clusters
tsne_clusters_interactive
#¬

#######################################################################

library(data.table)
library(griph)
library(scater)
library(AnnotationDbi)
library(stats)
library(org.Hs.eg.db)
library(scater)
library(ggplot2)
library(dplyr)
library(scran)
library(gridExtra)
library(tidyr)

# Load the filtered counts data
counts_filtered <- fread("re_analysis_counts_afterQC.csv")

# Remove rows with NaN and duplicates
counts_filtered <- na.omit(counts_filtered)  
counts_filtered <- unique(counts_filtered)

# Convert gene symbols to Entrez IDs
human <- org.Hs.eg.db
entrez_ids <- mapIds(human, keys = counts_filtered$external_gene_name, keytype = "SYMBOL", column = "ENTREZID")
counts_filtered$external_gene_name <- entrez_ids

# Extract the count matrix and set row names to Entrez IDs
count_matrix_filtered <- as.matrix(counts_filtered[, -1]) 
rownames(count_matrix_filtered) <- counts_filtered$external_gene_name

# Log-normalize the counts
log_counts_filtered <- log2(count_matrix_filtered + 1)

# Infer cell cycle stages
cc_scores <- predictCellCycle(log_counts_filtered, org = "human.Whitfield", cor_thr = 0.2, refine_iter = 200)

# Calculate library complexity (percentage of detected genes)
lib_complexity <- colSums(log_counts_filtered > 0) / nrow(log_counts_filtered) * 100

# Create a design matrix with cell cycle scores and library complexity
design_matrix <- model.matrix(~ cc_scores + lib_complexity)

# Regress out cell cycle and library complexity using lm
reg_counts <- t(apply(log_counts_filtered, 1, function(y) {
  fit <- lm(y ~ design_matrix)
  residuals(fit)
}))

# Add row and column names to the corrected counts matrix
rownames(reg_counts) <- rownames(log_counts_filtered)
colnames(reg_counts) <- colnames(log_counts_filtered)

# Update the SingleCellExperiment object with the new corrected counts
# First, subset the sce object to match the dimensions of reg_counts
sce_subset <- sce[rownames(sce) %in% rownames(reg_counts), colnames(sce) %in% colnames(reg_counts)]
assay(sce_subset, "logcounts") <- reg_counts[rownames(sce_subset), colnames(sce_subset)]
sce <- sce_subset

# Extract library complexity (percentage of detected genes)
sce$lib_complexity <- lib_complexity[colnames(sce)]

# Assign the cell names from sce to the cc_scores vector
names(cc_scores) <- colnames(sce)

# Subset cc_scores to match the column names of sce
cc_scores_subset <- cc_scores[colnames(sce)]

# Assign the subsetted cell cycle scores to the sce object
sce$cell_cycle_phase <- cc_scores_subset

# Categorize library complexity into four levels based on quartiles
lib_complexity_levels <- cut(sce$lib_complexity, breaks = quantile(sce$lib_complexity), 
                             labels = c("L1", "L2", "L3", "L4"), include.lowest = TRUE)
sce$lib_complexity_level <- lib_complexity_levels

# Plot t-SNE colored by library complexity levels (Fig. 2a equivalent)
png("tsne_libcomplexity_levels.png", width = 800, height = 600)
plotTSNE(sce, colour_by = "lib_complexity_level") + 
  ggtitle("t-SNE by Library Complexity Levels") +
  scale_color_manual(values = c("L1" = "blue", "L2" = "green", "L3" = "orange", "L4" = "red"))
dev.off()

#¬
# Create interactive t-SNE plot colored by library complexity levels
tsne_libcomplexity_interactive <- plotTSNE(sce, colour_by = "lib_complexity_level") +
  ggtitle("t-SNE by Library Complexity Levels") +
  scale_color_manual(values = c("L1" = "blue", "L2" = "green", "L3" = "orange", "L4" = "red"))

# Convert to interactive plotly object
tsne_libcomplexity_interactive <- ggplotly(tsne_libcomplexity_interactive, tooltip = "text")

# Adjust layout to include legend and title
tsne_libcomplexity_interactive <- tsne_libcomplexity_interactive %>%
  layout(title = "t-SNE by Library Complexity Levels",
         showlegend = TRUE,
         legend = list(x = 1.02, y = 0.5, orientation = "v"))

# Print the interactive plot
tsne_libcomplexity_interactive
#¬


# Combine cell cycle phases to match the original paper
sce$cell_cycle_phase_combined <- factor(sce$cell_cycle_phase, 
                                        levels = c("G1.S", "G1.S/S", "S", "S/G2", "G2", "G2/G2.M", "G2.M", "G2.M/M.G1", "M.G1", "M.G1/G1.S"),
                                        labels = c("G1/S", "G1/S", "S", "S/G2", "G2/M", "G2/M", "G2/M", "G1", "G1", "G1/S"))

# Plot t-SNE colored by combined cell cycle phases (Fig. 2b equivalent)
png("tsne_cellcycle_combined.png", width = 800, height = 600)
plotTSNE(sce, colour_by = "cell_cycle_phase_combined") + ggtitle("t-SNE by Combined Cell Cycle Phases")
dev.off()

#¬
# Combine cell cycle phases to match the original paper
sce$cell_cycle_phase_combined <- factor(sce$cell_cycle_phase,
                                        levels = c("G1.S", "G1.S/S", "S", "S/G2", "G2", "G2/G2.M", "G2.M", "G2.M/M.G1", "M.G1", "M.G1/G1.S"),
                                        labels = c("G1/S", "G1/S", "S", "S/G2", "G2/M", "G2/M", "G2/M", "G1", "G1", "G1/S"))


# Create a static t-SNE plot colored by combined cell cycle phases
tsne_cellcycle_combined <- plotTSNE(sce, colour_by = "cell_cycle_phase_combined") + 
  ggtitle("t-SNE by Combined Cell Cycle Phases")

# Convert the static plot to an interactive plotly object
tsne_cellcycle_combined_interactive <- ggplotly(tsne_cellcycle_combined, tooltip = "text")

# Print the interactive plot
tsne_cellcycle_combined_interactive
#¬




cell_cycle_distribution <- data.frame(
  Cluster = factor(sce$clusters_tsne),
  CellCycle = factor(sce$cell_cycle_phase_combined)
)

# Count the number of cells per cell cycle stage in each cluster
cell_cycle_distribution <- cell_cycle_distribution %>%
  group_by(Cluster, CellCycle) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  mutate(Percentage = Count / sum(Count) * 100)

# Create a data frame for library complexity distribution per cluster
library_complexity_distribution <- data.frame(
  Cluster = factor(sce$clusters_tsne),
  LibraryComplexity = factor(sce$lib_complexity_level)
)

# Count the number of cells per library complexity level in each cluster
library_complexity_distribution <- library_complexity_distribution %>%
  group_by(Cluster, LibraryComplexity) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  mutate(Percentage = Count / sum(Count) * 100)

# For the cell cycle bar graph
png("cell_cycle_bar_graph.png", width = 800, height = 600)
ggplot(cell_cycle_distribution, aes(x = Cluster, y = Percentage, fill = CellCycle)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  labs(x = "Cluster", y = "Cluster Composition (%)", fill = "Cell Cycle") +
  scale_fill_manual(values = c("G1" = "brown1", "G1/S" = "darkolivegreen3", "G2/M" = "cadetblue3", "S/G2" = "darkorchid2")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

# For the library complexity bar graph
png("library_complexity_bar_graph.png", width = 800, height = 600)
ggplot(library_complexity_distribution, aes(x = Cluster, y = Percentage, fill = LibraryComplexity)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  labs(x = "Cluster", y = "Cluster Composition (%)", fill = "Library Complexity") +
  scale_fill_manual(values = c("L1" = "brown1", "L2" = "darkolivegreen3", "L3" = "cadetblue3", "L4" = "darkorchid2")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()

#¬
# For the cell cycle bar graph
cell_cycle_bar_graph_interactive <- ggplot(cell_cycle_distribution, aes(x = Cluster, y = Percentage, fill = CellCycle, text = paste("Cluster:", Cluster, "<br>Cell Cycle:", CellCycle, "<br>Percentage:", round(Percentage, 2), "%"))) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  labs(x = "Cluster", y = "Cluster Composition (%)", fill = "Cell Cycle") +
  scale_fill_manual(values = c("G1" = "brown1", "G1/S" = "darkolivegreen3", "G2/M" = "cadetblue3", "S/G2" = "darkorchid2")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

cell_cycle_bar_graph_interactive <- ggplotly(cell_cycle_bar_graph_interactive, tooltip = "text")
cell_cycle_bar_graph_interactive

# For the library complexity bar graph
library_complexity_bar_graph_interactive <- ggplot(library_complexity_distribution, aes(x = Cluster, y = Percentage, fill = LibraryComplexity, text = paste("Cluster:", Cluster, "<br>Library Complexity:", LibraryComplexity, "<br>Percentage:", round(Percentage, 2), "%"))) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  labs(x = "Cluster", y = "Cluster Composition (%)", fill = "Library Complexity") +
  scale_fill_manual(values = c("L1" = "brown1", "L2" = "darkolivegreen3", "L3" = "cadetblue3", "L4" = "darkorchid2")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

library_complexity_bar_graph_interactive <- ggplotly(library_complexity_bar_graph_interactive, tooltip = "text")
library_complexity_bar_graph_interactive
# Print the interactive plots
cell_cycle_bar_graph_interactive
library_complexity_bar_graph_interactive
#¬




# Create a data frame for model distribution per cluster
model_distribution <- data.frame(
  Cluster = factor(sce$clusters_tsne),
  Model = factor(sce$model)
)

# Count the number of cells per model in each cluster
model_distribution <- model_distribution %>%
  group_by(Cluster, Model) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  mutate(Percentage = Count / sum(Count) * 100)

# Plot the model distribution bar graph
png("model_distribution_bar_graph.png", width = 800, height = 600)
ggplot(model_distribution, aes(x = Cluster, y = Percentage, fill = Model)) +
  geom_bar(stat = "identity", position = "fill") +
  theme_minimal() +
  labs(x = "Cluster", y = "Cluster Composition (%)", fill = "Model") +
  scale_fill_manual(values = c("PDX1" = "brown1", "PDX2" = "darkolivegreen3", "PDX3" = "cadetblue3", "PDX4" = "darkorchid2", "MDAMB231" = "pink")) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))
dev.off()
