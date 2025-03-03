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

# Plot results (MA-plot)
plotMA(res, ylim = c(-2, 2))

# Save the plot
dev.copy(png, “/Users/alannahodonnell/Downloads/salmon_out_M/MA_plot_.png”)
dev.off()
library(ggplot2)
# PCA
vsd <- vst(dds, blind = FALSE)  # Variance stabilizing transformation
pca_data <- plotPCA(vsd, intgroup = “condition”, returnData = TRUE)

# Print PCA data
print(head(pca_data))
#Save PCA plot
pca_plot <- ggplot(pca_data, aes(x = PC1, y = PC2, color = condition)) +
  geom_point(size = 3) +
  xlab(paste0(“PC1: “, round(attr(pca_data, “percentVar”)[1], 1), “% variance”)) +
  ylab(paste0(“PC2: “, round(attr(pca_data, “percentVar”)[2], 1), “% variance”)) +
  theme_minimal() +
  ggtitle(“PCA of RNA-seq Data”)

# Save the PCA plot as a PNG file
ggsave(“/Users/alannahodonnell/Downloads/salmon_out/PCA_plot.png”, plot = pca_plot)

# Display the PCA plot
print(pca_plot)
