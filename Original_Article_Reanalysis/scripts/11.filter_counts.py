mport pandas as pd

# Read the featureCounts output
counts_df = pd.read_csv('counts_updated_with_SampleTitles.csv', index_col=0)
print("Original dataframe shape:", counts_df.shape)
print("Example data:\n", counts_df.iloc[:5, :5])

# Step 1: Calculate the total number of reads mapping to the human transcriptome for each library
total_reads_per_lib = counts_df.sum()
print("\nTotal reads per library (first 5):\n", total_reads_per_lib.head())

# Step 2: Filter out libraries with <= 100k reads mapping to the human transcriptome
filtered_df = counts_df.loc[:, total_reads_per_lib > 100000]
print("\nDataframe shape after filtering libraries with <= 100k reads:", filtered_df.shape)
print("Example data after filtering libraries with <= 100k reads:\n", filtered_df.iloc[:5, :5])

# Step 3: Calculate the number of expressed genes for each remaining library
num_expressed_genes_per_lib = (filtered_df > 0).sum()
print("\nNumber of expressed genes per library (first 5):\n", num_expressed_genes_per_lib.head())

# Step 4: Filter out libraries expressing less than 2,346 genes or more than 9,884 genes
min_genes = 2346
max_genes = 9884
filtered_df = filtered_df.loc[:, (num_expressed_genes_per_lib >= min_genes) & (num_expressed_genes_per_lib <= max_genes)]
print("\nDataframe shape after filtering libraries expressing <2,346 or >9,884 genes:", filtered_df.shape)
print("Example data after filtering libraries expressing <2,346 or >9,884 genes:\n", filtered_df.iloc[:5, :5])

# Step 5: Check if the remaining cells express 10%-40% of the total number of unique genes observed (23,459)
total_genes_observed = 23459
min_percent = 0.1
max_percent = 0.4

num_cells_expr_10_40_percent = ((filtered_df > 0).sum() >= min_percent * total_genes_observed) & ((filtered_df > 0).sum() <= max_percent * total_genes_observed)
print(f"\nNumber of cells expressing 10%-40% of total unique genes observed: {num_cells_expr_10_40_percent.sum()}")
print(f"Percentage of cells expressing 10%-40% of total unique genes observed: {num_cells_expr_10_40_percent.mean():.2%}")

print("\nFinal dataframe shape:", filtered_df.shape)
print("Example data from final dataframe:\n", filtered_df.iloc[:5, :5])

# Save the filtered counts matrix
filtered_df.to_csv('filtered_counts.csv')
print("\nFiltered counts matrix saved to 'filtered_counts.csv'")
