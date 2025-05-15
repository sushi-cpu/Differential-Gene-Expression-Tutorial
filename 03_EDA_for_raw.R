library(affy)
library(hgu133plus2.db)
library(AnnotationDbi)
library(pheatmap)
library(ggplot2)

# 01_Load CEL files
cel_files <- list.celfiles("data/raw/", full.names = TRUE) # your path
raw_data <- ReadAffy(filenames = cel_files)

# 02_RMA Normalization
norm_data <- rma(raw_data)

# 03_Extract Expression Matrix
expr_matrix <- exprs(norm_data)

# 04_Define sample groups manually
sample_labels <- factor(c(rep("Tumor", 54), rep("Normal", 12)))
if(length(sample_labels) != ncol(expr_matrix)) {
  stop("Length of sample_labels does not match number of samples!")
}

# 05_PCA on normalized data
pca <- prcomp(t(expr_matrix), scale. = TRUE)
pca_data <- data.frame(pca$x)
pca_data$Group <- sample_labels

ggplot(pca_data, aes(x = PC1, y = PC2, color = Group)) +
  geom_point(size = 3) +
  theme_minimal() +
  labs(title = "PCA from CEL files", color = "Sample Type")

# 06_Top 50 most variable genes (adjust number here)
top_var_genes <- order(apply(expr_matrix, 1, var), decreasing = TRUE)[1:10]
top_genes <- expr_matrix[top_var_genes, ]
top_genes <- as.data.frame(top_genes)

# 07_Map probe IDs to gene symbols
probe_ids <- rownames(top_genes)
gene_symbols <- mapIds(hgu133plus2.db,
                       keys = probe_ids,
                       column = "SYMBOL",
                       keytype = "PROBEID",
                       multiVals = "first")

# 08_Annotate and clean
top_genes$GeneSymbol <- gene_symbols
top_genes <- top_genes[!is.na(top_genes$GeneSymbol), ]
top_genes <- top_genes[!duplicated(top_genes$GeneSymbol), ]
rownames(top_genes) <- top_genes$GeneSymbol
top_genes$GeneSymbol <- NULL

# 09_Annotation for heatmap columns
annotation_col <- data.frame(Group = sample_labels)
rownames(annotation_col) <- colnames(top_genes)

# 10_Plot heatmap
pheatmap(top_genes,
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         annotation_col = annotation_col,
         main = "Top 10 Variable Genes (From CEL Files)",
         filename = "results/Top 10 Variable Genes (From CEL Files).png",
         width = 10,      # width in inches
         height = 6,     # height in inches
         units = "in",   # units for width and height: "in", "cm", or "px"
         dpi = 600)      # resolution for raster formats like PNG)


