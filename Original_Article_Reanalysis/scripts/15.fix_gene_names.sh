# Step 1: Extract the first column from re_analysis_counts_afterQC.csv
cut -d, -f1 re_analysis_counts_afterQC.csv > first_column.csv

# Step 2: Remove the first column from re_analysis_CCcorrected_logcounts.csv and keep the rest
cut -d, -f2- re_analysis_CCcorrected_logcounts.csv > remaining_columns.csv

# Step 3: Combine the new first column with the remaining columns and overwrite the original file
paste -d, first_column.csv remaining_columns.csv > re_analysis_CCcorrected_logcounts.csv

# Step 4: Clean up the temporary files
rm first_column.csv remaining_columns.csv
