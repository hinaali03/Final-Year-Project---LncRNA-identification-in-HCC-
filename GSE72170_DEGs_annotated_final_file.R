# Step 1: Load the necessary data files
combined_data <- read.csv("GSE72170_combined_data_with_SEQ.csv", stringsAsFactors = FALSE)
annotation_data <- read.csv("GSE72170_lncRNA_annotation_file.csv", stringsAsFactors = FALSE)

# Step 2: Merge the data frames based on the 'ID' column
merged_data <- merge(combined_data, annotation_data[, c("ID", "lncRNA_gene_id", "lncRNA_gene_name")], by = "ID", all.x = TRUE)

# Step 3: Remove rows where 'lncRNA_gene_id' or 'lncRNA_gene_name' are NA
merged_data <- merged_data[!is.na(merged_data$lncRNA_gene_id) & !is.na(merged_data$lncRNA_gene_name), ]

# Step 4: Check the merged data
head(merged_data)

# Optionally, save the cleaned merged data to a new CSV file
write.csv(merged_data, "Merged_GSE72170_cleaned_data.csv", row.names = FALSE)

# Step 4: Subset the final merged data to include only 'lncRNA_gene_id', 'lncRNA_gene_name', and 'AveExpr'
final_data <- merged_data[, c("lncRNA_gene_id", "lncRNA_gene_name", "AveExpr")]

# Step 5: Save the final data to a new CSV file
write.csv(final_data, "GSE72170_final_DEGs.csv", row.names = FALSE)
