awk -F, '!($1 ~ /NA_/ || $1 ~ /_/) {print}' 1.csv > 2.csv
sed '1s/^V1,/external_gene_name,/' 2.csv > re_analysis_counts_afterQC.csv
