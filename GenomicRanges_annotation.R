# Install necessary packages if not already installed
install.packages("GenomicRanges")
install.packages("rtracklayer")
install.packages("dplyr")



# Load the necessary libraries
library(dplyr)
library(GenomicRanges)
library(rtracklayer)
library(dplyr)

## Making the sequence fasta file

# 1. Load the DEGs file (adjust the file path and column names accordingly)
degs_file <- "GSE72170_combined_data_with_SEQ.csv"
degs <- read.csv(degs_file)

# 2. Extract the first 10 sequences and their corresponding IDs
sequences <- degs$SEQUENCE
ids <- degs$ID  # Adjust this to your actual ID column name

# 3. Save these sequences to a FASTA file for BLAST
fasta_file <- "deg_sequences.fasta"

write_fasta <- function(sequences, ids, file_name) {
  sink(file_name)
  for (i in 1:length(sequences)) {
    cat(paste0(">", ids[i], "\n", sequences[i], "\n"))
  }
  sink()
}

write_fasta(sequences, ids, fasta_file)

## Blast on windows subsystem linux
#command:
# blastn -query deg_sequences.fasta -db GRCh38_p14_db -outfmt 6 -out blast_results.txt -num_threads 4



## Annotating Blast result

# 1. Read the BLAST results into R
blast_results <- read.table("blast_results.txt", header = FALSE)

# Assign appropriate column names to the BLAST results
colnames(blast_results) <- c("query_id", "subject_id", "percent_identity", "alignment_length", 
                             "mismatches", "gap_openings", "q_start", "q_end", "s_start", "s_end", 
                             "e_value", "bit_score")

# 2. Read the GFF3 file for long non-coding RNA annotations
lncRNA_gff <- import("gencode.v47.annotation.gff3")



# Ensure start is always less than or equal to end for each BLAST hit
blast_results$start_fixed <- pmin(blast_results$s_start, blast_results$s_end)
blast_results$end_fixed <- pmax(blast_results$s_start, blast_results$s_end)

# Create GenomicRanges for the BLAST hits with fixed start and end positions
blast_granges <- GRanges(
  seqnames = blast_results$subject_id,  # subject_id corresponds to the genome chromosome
  ranges = IRanges(start = blast_results$start_fixed, end = blast_results$end_fixed),
  query_id = blast_results$query_id
)

# Inspect the created GRanges object
head(blast_granges)

######################################################
# 4. Find Overlaps Between BLAST Hits and lncRNA Annotations
overlaps <- findOverlaps(blast_granges, lncRNA_gff)

# View the overlaps
head(overlaps)

# 5. Extract the Overlapping Annotations (lncRNAs) from the GFF3 File
lncRNA_hits <- lncRNA_gff[subjectHits(overlaps)]



# View the annotated BLAST results
head(blast_results)


############ debugging:##############################################################3
# Check the number of overlaps
length(subjectHits(overlaps))  # This should give us the number of valid overlaps
# Check the first few overlaps
head(subjectHits(overlaps))
head(queryHits(overlaps))
###########################



# 6. Check the metadata structure for gene_id and gene_name
#metadata <- elementMetadata(lncRNA_hits)
#head(metadata)  # This will show the first few rows and columns of metadata



# Initialize columns with NA values in case there is no overlap
blast_results$lncRNA_gene_id <- NA
blast_results$lncRNA_gene_name <- NA

# Check the overlaps between blast_results and lncRNA_hits
for (i in 1:nrow(blast_results)) {
  # Find the corresponding overlaps for each BLAST hit (query)
  hit_idx <- which(queryHits(overlaps) == i)
  
  if (length(hit_idx) > 0) {
    # Check if multiple overlaps exist and select one or handle accordingly
    blast_results$lncRNA_gene_id[i] <- elementMetadata(lncRNA_hits)$gene_id[hit_idx[1]]  # Selecting the first overlap
    blast_results$lncRNA_gene_name[i] <- elementMetadata(lncRNA_hits)$gene_name[hit_idx[1]]  # Selecting the first overlap
    blast_results$lncRNA_gene_type[i] <- elementMetadata(lncRNA_hits)$gene_type[hit_idx[1]]  # Selecting the first overlap
  }
}

# View the annotated BLAST results (you can check the first few rows)
head(blast_results)

#Sort by query_id and alignment_length (in descending order) to keep the best (longest) hit for each query_id
best_hits <- blast_results %>%
  arrange(query_id, desc(alignment_length)) %>%
  distinct(query_id, .keep_all = TRUE)

# View the first few rows of the cleaned-up results
head(best_hits)

# Save the final annotated results to a CSV file
write.csv(best_hits, "best_annotated_blast_results.csv", row.names = FALSE)

# Remove rows where 'lncRNA_gene_name' is NA
best_hits_no_na <- best_hits[!is.na(best_hits$lncRNA_gene_name), ]


# Save the final annotated results to a CSV file without NAs in the 'lncRNA_gene_name' column
write.csv(best_hits_no_na, "best_annotated_blast_results_noNA.csv", row.names = FALSE)

# 7. Save the Annotated Results to a CSV File
write.csv(blast_results, "annotated_blast_results.csv", row.names = FALSE)
