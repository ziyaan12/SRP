# Set working directory (Session > Set Working Directory > To Source File Location)

# Load the filtered counts data
counts <- fread("filtered_counts.csv")

# Remove rows with NaN values
counts <- counts[complete.cases(counts), ]

# Remove duplicate rows
counts <- unique(counts)

# Save the cleaned data as a new CSV file
fwrite(counts, "1.csv")
