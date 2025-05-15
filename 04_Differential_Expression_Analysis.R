# Load required libraries
library(affy)
library(limma)
library(hgu133plus2.db)
library(ggplot2)
library(AnnotationDbi)
library(pheatmap)
library(ggrepel)
library(plotly)
library(heatmaply)
library(htmlwidgets)

# -----------------------------------
# Step 1: Read and Normalize the CEL Files
# -----------------------------------
cel_files <- list.celfiles("data/raw/", full.names = TRUE)
affy_data <- ReadAffy(filenames = cel_files)
eset <- rma(affy_data)
exprs_matrix <- exprs(eset)

# -----------------------------------
# Step 2: Define Experimental Groups
# -----------------------------------
group <- factor(c(rep("HER2_Pos", 28), rep("HER2_Neg", 26), rep("Normal", 12)))
design <- model.matrix(~0 + group)
colnames(design) <- levels(group)

contrast.matrix <- makeContrasts(
  HER2_Pos_vs_Neg = HER2_Pos - HER2_Neg,
  HER2_Pos_vs_Normal = HER2_Pos - Normal,
  HER2_Neg_vs_Normal = HER2_Neg - Normal,
  levels = design
)

# -----------------------------------
# Step 3: Fit Linear Model and Apply Contrasts
# -----------------------------------
fit <- lmFit(exprs_matrix, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)

# -----------------------------------
# Step 4: Extract DEGs for HER2_Pos vs Normal
# -----------------------------------
deg_results <- topTable(fit2, coef = "HER2_Pos_vs_Normal", adjust.method = "fdr", number = Inf)

# -----------------------------------
# Step 5: Annotate with Gene Symbols
# -----------------------------------
deg_results$GeneSymbol <- mapIds(hgu133plus2.db,
                                 keys = rownames(deg_results),
                                 column = "SYMBOL",
                                 keytype = "PROBEID",
                                 multiVals = "first")

# Clean up and deduplicate
deg_results_annot <- na.omit(deg_results)
deg_results_annot_unique <- deg_results_annot[!duplicated(deg_results_annot$GeneSymbol), ]
deg_results_annot_unique$ProbeID <- rownames(deg_results_annot_unique)
rownames(deg_results_annot_unique) <- deg_results_annot_unique$GeneSymbol
deg_results_annot_unique$GeneSymbol <- NULL

# -----------------------------------
# Step 6: Add DEG Status and Save Results
# -----------------------------------
# Define thresholds
logFC_cutoff <- 1
pval_cutoff <- 0.05

# Create DEG status labels
deg_results_annot_unique$gene <- rownames(deg_results_annot_unique)
deg_results_annot_unique$status <- with(deg_results_annot_unique, ifelse(
  adj.P.Val < pval_cutoff & logFC > logFC_cutoff, "Upregulated",
  ifelse(adj.P.Val < pval_cutoff & logFC < -logFC_cutoff, "Downregulated", "Not Significant")
))

# Reorder columns
deg_results_annot_unique <- deg_results_annot_unique[, c("ProbeID", 
                                                         "logFC", 
                                                         "AveExpr", 
                                                         "t", 
                                                         "P.Value", 
                                                         "adj.P.Val", 
                                                         "B", 
                                                         "gene", 
                                                         "status")]

# Save to CSV
write.csv(deg_results_annot_unique, "DEGs.csv", row.names = TRUE)

# -----------------------------------
# Step 7: Interactive Volcano Plot
# -----------------------------------
status_colors <- c("Upregulated" = "firebrick", 
                   "Downregulated" = "royalblue", 
                   "Not Significant" = "grey")

p <- plot_ly(deg_results_annot_unique,
        x = ~logFC,
        y = ~-log10(adj.P.Val),
        text = ~gene,
        color = ~status,
        colors = status_colors,
        mode = "markers",
        type = "scatter",
        marker = list(size = 5)) %>%
  layout(title = "Interactive Volcano Plot: HER2_Pos vs Normal",
         xaxis = list(title = "Log2 Fold Change"),
         yaxis = list(title = "-log10 Adjusted P-Value"))

saveWidget(p, "results/volcano_plot.html", selfcontained = TRUE)

# -----------------------------------
# Step 8: Heatmap of Top 10 DEGs
# -----------------------------------
# Subset top significant genes
significant_genes <- subset(deg_results_annot_unique, 
                            adj.P.Val < pval_cutoff & abs(logFC) > logFC_cutoff)
top_genes <- significant_genes[order(-abs(significant_genes$logFC)), ]
top_genes <- head(top_genes, 10)

# Get expression data
exprs_sig <- exprs_matrix[top_genes$ProbeID, ]
rownames(exprs_sig) <- rownames(top_genes)

# Sample annotation
sample_annotation <- data.frame(Group = group)
rownames(sample_annotation) <- colnames(exprs_sig)

# Interactive heatmap
s <- heatmaply(exprs_sig,
          colors = viridis::viridis(100),
          Rowv = TRUE,
          Colv = TRUE,
          col_side_colors = sample_annotation,
          scale = "row",
          main = "Top 10 Most Significant DEGs",
          showticklabels = c(TRUE, FALSE),
          fontsize_row = 8,
          dendrogram = "both",
          labCol = rep("", ncol(exprs_sig)))

saveWidget(s, "results/heatmaply_top10_DEGs.html", selfcontained = TRUE)

# Static heatmap
pheatmap(exprs_sig,
         annotation_col = sample_annotation,
         scale = "row",
         show_rownames = TRUE,
         show_colnames = FALSE,
         fontsize_row = 8,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         main = "Heatmap of Top 10 DEGs",
         filename = "results/Heatmap of Top 10 DEGs.png",
         width = 10,      # width in inches
         height = 6,     # height in inches
         units = "in",   # units for width and height: "in", "cm", or "px"
         dpi = 600)      # resolution for raster formats like PNG)

