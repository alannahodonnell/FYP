if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("DESeq2")

BiocManager::install("tximport")

install.packages("ggplot2")



#Define the directory path
dir_path <- “/Users/alannahodonnell/Downloads/salmon_out/”

# List all subdirectories
subdirs <- list.dirs(dir_path, full.names = TRUE, recursive = FALSE)

# List all quant.sf files within each subdirectory
quant_files <- unlist(lapply(subdirs, function(subdir) {
  list.files(subdir, pattern = “quant.sf”, full.names = TRUE)
}))

# Print the files to check if they are correctly listed
print(quant_files)
# Check if the number of files is as expected
if (length(quant_files) != 8) {
  stop(“The number of quant.sf files found does not match the expected number (8). Check the directory and pattern.”)
}

sample_names <- c(“SRR6175537”, “SRR6175538”, “SRR6175539”, “SRR6175540”, “SRR6175541”, “SRR6175542”, “SRR6175543”, “SRR6175544”)

# Create a sample table
sample_table <- data.frame(
  sampleName = sample_names,
   ilename = quant_files,
  condition = factor(rep(c(“control”, “ko”), each = 4))
)

# Print the sample table to ensure it’s correct
print(sample_table)

# Read in Salmon quantification data
txi <- tximport(quant_files, type = “salmon”, txOut = TRUE)

# Create DESeq2 dataset
dds <- DESeqDataSetFromTximport(txi,
                                colData = sample_table,
                                design = ~ condition)

# Pre-filtering (optional but recommended)
dds <- dds[ rowSums(counts(dds)) > 1, ]

# Run the DESeq2 pipeline
dds <- DESeq(dds)

# Results
res <- results(dds)

# Output results to a file
write.csv(as.data.frame(res), 
          file = “/Users/alannahodonnell/Downloads/salmon_out/deseq2_results_M.csv”)

# Load necessary libraries
library(ggplot2)
library(ggrepel)

# Read DESeq2 results (adjust file path if needed)
res_df <- read.csv("/Users/alannahodonnell/Downloads/.csv", sep = ",", stringsAsFactors = FALSE)

# Fix column names if they contain dots instead of underscores
colnames(res_df) <- gsub("\\.", "_", colnames(res_df))

# Remove rows where pvalue or padj are NA
res_df <- res_df[!is.na(res_df$pvalue) & !is.na(res_df$padj), ]

# Ensure necessary columns are numeric
res_df$pvalue <- as.numeric(res_df$pvalue)
res_df$padj <- as.numeric(res_df$padj)
res_df$log2FoldChange <- as.numeric(res_df$log2FoldChange)

# Calculate -log10(p-value) for better visualization
res_df$logP <- -log10(res_df$pvalue)

# Define significance categories
res_df$significance <- ifelse(
  res_df$padj < 0.05 & res_df$log2FoldChange > 1, "Upregulated",
  ifelse(res_df$padj < 0.05 & res_df$log2FoldChange < -1, "Downregulated", "Not Significant")
)

# Select top 10 upregulated and downregulated genes for labeling
upregulated_genes <- res_df[res_df$significance == "Upregulated", ]
downregulated_genes <- res_df[res_df$significance == "Downregulated", ]

# Order by padj and select top 10 for each
top_upregulated_genes <- upregulated_genes[order(upregulated_genes$padj), ]
top_downregulated_genes <- downregulated_genes[order(downregulated_genes$padj), ]

top_upregulated_genes <- head(top_upregulated_genes, 10)
top_downregulated_genes <- head(top_downregulated_genes, 10)

# Combine top upregulated and downregulated genes
top_genes <- rbind(top_upregulated_genes, top_downregulated_genes)

# Create Volcano Plot
volcano_plot <- ggplot(res_df, aes(x = log2FoldChange, y = logP, color = significance)) +
  geom_point(alpha = 0.7) +
  scale_color_manual(values = c("Upregulated" = "red", "Downregulated" = "blue", "Not Significant" = "gray")) +
  theme_minimal() +
  labs(title = "Volcano Plot", x = "Log2 Fold Change", y = "-Log10 P-value") +
  theme(legend.title = element_blank()) +
  geom_text_repel(data = top_genes, aes(label = Gene_Symbol),  # Ensure correct column name
                  size = 5, box.padding = 0.5, max.overlaps = 15) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +  # Add threshold line
  xlim(c(-10, 10)) +  # Set x-axis range from -10 to 10
  guides(color = guide_legend(override.aes = list(size = 4)))  # Adjust legend without letters

# Save the plot as PNG
png("/Users/alannahodonnell/Downloads/Volcano_plot_8.png", width = 8, height = 6, units = "in", res = 300)
print(volcano_plot)
dev.off()
