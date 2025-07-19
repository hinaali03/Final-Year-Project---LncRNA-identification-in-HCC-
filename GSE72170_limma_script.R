# Install GEOquery package if you haven't already
#install.packages("BiocManager")
#BiocManager::install("GEOquery")

# Load the GEOquery package
library(GEOquery)
library(limma)


# Step 1: Download the dataset (Example using GEO database)
# Replace 'GSEXXXX' with your dataset's GEO accession number
geo_accession <- "GSE72170"  # Replace this with your GEO ID

gse <- getGEO(geo_accession, GSEMatrix = TRUE)


# Step 2: Extract the expression data matrix (assuming it's in the first element)
expr_data <- exprs(gse[[1]])

# Step 3: Create a phenotype data frame
# This assumes the samples are labeled "Tumor" and "Normal"
pheno_data <- pData(gse[[1]])
head(pheno_data)
#tail(pheno_data)
#colnames(pheno_data)
#Save the phenotype data (pheno_data) to a CSV file
write.csv(pheno_data, file = "pheno_data.csv", row.names = FALSE)

# Step 4: Check the tissue column
# Let's confirm the unique values in the tissue column
table(pheno_data$`tissue:ch1`)


# Step 8: Modify the column to create a factor for Tumor vs Noncancerous
pheno_data$sample_type <- factor(pheno_data$`tissue:ch1`, levels = c("Tumor tissue", "Matched noncancerous tissue"),
                                 labels = c("Tumor", "Normal"))

# Verify the modified class of the samples (Tumor vs Normal)
table(pheno_data$sample_type)

# Step 9: Normalize the expression data (optional, but recommended)
expr_data_normalized <- normalizeBetweenArrays(expr_data)

# Step 10: Create a design matrix for the linear model
design <- model.matrix(~ 0 + factor(pheno_data$sample_type))
colnames(design) <- levels(factor(pheno_data$sample_type))

# Step 11: Fit the linear model using limma
fit <- lmFit(expr_data_normalized, design)

# Step 12: Set up contrasts to compare Tumor vs Normal
contrast.matrix <- makeContrasts(Tumor - Normal, levels = design)

# Step 13: Apply the contrasts and calculate statistics
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# Step 14: View the top differentially expressed genes (DEGs)
results <- topTable(fit2, adjust = "fdr", number = Inf)  # Adjust for multiple testing (FDR)
head(results)

# Apply p-value threshold of 0.05 to filter DEGs
filtered_results <- results[results$P.Value < 0.05, ]

# Step 15: Save filtered results to CSV
write.csv(filtered_results, file = "GSE72170_DEGs_results_p05.csv")

# Step 16: Save the complete results to CSV
write.csv(results, file = "DEGs_results.csv")

# Plot the volcano plot to visualize DEGs
volcanoplot(fit2, coef = 1, highlight = 10, main = "Volcano Plot for Tumor vs Normal")

# Optional: Save the volcano plot as a PNG file
png("volcano_plot.png", width = 800, height = 600)
volcanoplot(fit2, coef = 1, highlight = 10, main = "Volcano Plot for Tumor vs Normal")
dev.off()

########### Annotation with GB_ACC

DEGs <- read.csv("GSE72170_DEGs_results_p05.csv")
head(DEGs)


# Step 2: Use the DEGs data as the expression data
# Assuming DEGs has a column for identifiers (e.g., "ID") and expression data for each sample

# Rename the first column as 'ID' if necessary, and ensure the row names are in the first column
deg_expr <- DEGs
#expr_data <- cbind(ID = rownames(expr_data), expr_data)
head(deg_expr)

# Save the first few rows to a CSV file
#write.csv(head(expr_data, 10), "expr_data_with_ID.csv", row.names = FALSE)

# Convert expr_data to a data frame (if it's not already)
deg_expr <- as.data.frame(deg_expr)

# Step 3: Load the SOFT formatted family file (GSE84004_family.soft) that you've downloaded
soft_file_path <- "GSE72170_family.soft"  # Specify the path to your downloaded SOFT file
soft_data <- read.table(soft_file_path, sep = "\t", header = TRUE, quote = "\"", comment.char = "!")

# Assuming the feature column is named 'Feature' (adjust the column name as needed)
# Filter the data to keep only rows where the feature column contains 'ncRNA'
#soft_data <- soft_data[soft_data$Feature == "ncRNA", ]

# Step 5: Extract GB_ACC and corresponding identifiers (assuming it's in a column called "GB_ACC")
soft_data <- soft_data[, c("ID", "Sequence")]
head(soft_data)
print(nrow(soft_data))
write.csv(soft_data, "soft_data.csv", row.names = FALSE)



## GETTING SEQ
soft_data_SEQ <- soft_data[!is.na(soft_data$Sequence) & soft_data$Sequence != "", ]
print(nrow(soft_data_SEQ))
write.csv(soft_data, "soft_data_SEQ.csv", row.names = FALSE)

# Step 6: Perform the merge using dplyr's left_join (ensure matching IDs)
combined_expr_data <- merge(soft_data_SEQ, deg_expr, by = "ID")
sum(is.na(combined_expr_data$Sequence) | combined_expr_data$Sequence == "")
print(nrow(combined_expr_data))
# Step 7: Save the combined result to CSV (this will include GB_ACC)
write.csv(combined_expr_data, "combined_data_with_SEQ.csv", row.names = FALSE)
