# Load required libraries
library(GEOquery)
library(ggplot2)
library(pheatmap)
library(hgu133plus2.db)
library(AnnotationDbi)

# -----------------------------------
# Step 1: Read Processed Data
# -----------------------------------
processed_data <- read.table("data/processed/GSE29431_series_matrix.txt", 
                             header = TRUE, 
                             sep = "\t", 
                             row.names = 1, 
                             fill = TRUE, 
                             comment.char = "!", 
                             check.names = FALSE)

# -----------------------------------
# Step 2: Define Sample Groups Manually
# (Adjust numbers according to your dataset)
# -----------------------------------
group_labels <- factor(c(rep("Tumor", 54), rep("Normal", 12)))

# Check if group length matches number of samples
if(length(group_labels) != ncol(processed_data)){
  stop("Error: Length of group_labels does not match number of samples in processed_data!")
}

# -----------------------------------
# Step 3: Perform PCA on Samples
# -----------------------------------
pca <- prcomp(t(processed_data), scale. = TRUE)
pca_data <- data.frame(pca$x)
pca_data$Group <- group_labels

# PCA plot
ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  labs(title = "PCA of Processed Data", color = "Group") +
  theme_minimal()

# -----------------------------------
# Step 4: Select Top 10 Most Variable Genes
# -----------------------------------
top_var_genes <- order(apply(processed_data, 1, var), decreasing = TRUE)[1:10]
top_genes <- processed_data[top_var_genes, ]

# -----------------------------------
# Step 5: Map Probe IDs to Gene Symbols
# -----------------------------------
probe_ids <- rownames(top_genes)
gene_symbols <- mapIds(hgu133plus2.db,
                       keys = probe_ids,
                       column = "SYMBOL",
                       keytype = "PROBEID",
                       multiVals = "first")

# Add gene symbols to data frame
top_genes_annotated <- top_genes
top_genes_annotated$GeneSymbol <- gene_symbols

# Remove rows with missing gene symbols
top_genes_annotated <- top_genes_annotated[!is.na(top_genes_annotated$GeneSymbol), ]

# Remove duplicated gene symbols, keep first
top_genes_annotated <- top_genes_annotated[!duplicated(top_genes_annotated$GeneSymbol), ]

# Move gene symbols to rownames and remove the column
rownames(top_genes_annotated) <- top_genes_annotated$GeneSymbol
top_genes_annotated$GeneSymbol <- NULL

# -----------------------------------
# Step 6: Create Sample Annotation for Heatmap
# -----------------------------------
annotation_col <- data.frame(Group = group_labels)
rownames(annotation_col) <- colnames(processed_data)

# -----------------------------------
# Step 7: Plot Heatmap of Top Variable Genes
# -----------------------------------
pheatmap(top_genes_annotated, 
         cluster_rows = TRUE, 
         cluster_cols = TRUE, 
         annotation_col = annotation_col,
         show_rownames = TRUE, 
         main = "Heatmap of Top 10 Variable Genes",
         filename = "results/top_variable_genes_heatmap.png",
         width = 10,      # width in inches
         height = 6,     # height in inches
         units = "in",   # units for width and height: "in", "cm", or "px"
         dpi = 600)      # resolution for raster formats like PNG)


