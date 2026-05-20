# ============================================================
# Microarray Differential Expression Analysis — GSE233849
# Platform : Clariom S Human
# Comparison: MI vs Vehicle (5 vs 5 samples)
# Author   : [Your Name]
# Date     : 2025
# ============================================================
# This script performs end-to-end analysis of Affymetrix
# Clariom S Human microarray data including:
#   1. Quality control & normalization (RMA)
#   2. Differential expression analysis (limma)
#   3. Visualization (PCA, heatmap, volcano plot)
#   4. Gene network construction (WGCNA / igraph)
# ============================================================


# ============================================================
# SECTION 1 — INSTALL PACKAGES (run once)
# ============================================================

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install(c(
  "affy", "gcrma", "limma", "oligo",
  "pd.clariom.s.human",
  "clariomshumantranscriptcluster.db",
  "AnnotationDbi",
  "WGCNA"
))

install.packages(c("pheatmap", "igraph", "corrplot", "ggplot2", "ggrepel"))


# ============================================================
# SECTION 2 — LOAD LIBRARIES
# ============================================================

library(affy)
library(gcrma)
library(limma)
library(oligo)
library(pd.clariom.s.human)
library(clariomshumantranscriptcluster.db)
library(AnnotationDbi)
library(pheatmap)
library(ggplot2)
library(ggrepel)
library(WGCNA)
library(igraph)
library(corrplot)


# ============================================================
# SECTION 3 — SET WORKING DIRECTORY & READ CEL FILES
# ============================================================

setwd("GSE233849_RAW")   # folder containing all .CEL files
list.files()

# Read all CEL files in the current directory
data <- read.celfiles(list.celfiles())
data


# ============================================================
# SECTION 4 — QUALITY CONTROL (BEFORE NORMALIZATION)
# ============================================================

# Box-plot of raw probe intensities — samples should align after RMA
boxplot(data,
        main = "Raw Data — Probe Intensity Distribution",
        ylab = "log2 Intensity",
        las  = 2,
        col  = "lightblue")


# ============================================================
# SECTION 5 — NORMALIZATION (RMA)
# ============================================================
# RMA (Robust Multichip Average) performs background correction,
# quantile normalization, and log2 summarization at the probeset level.

norm_data <- rma(data)
norm_data

# Expression matrix: rows = probesets, columns = samples
expr_matrix <- exprs(norm_data)
dim(expr_matrix)
head(expr_matrix)

# Box-plot after normalization — distributions should now overlap
boxplot(expr_matrix,
        main = "Normalized Data — RMA",
        ylab = "log2 Intensity",
        las  = 2,
        col  = "lightgreen")


# ============================================================
# SECTION 6 — ANNOTATE PROBE IDs → GENE SYMBOLS
# ============================================================
# Map Affymetrix probeset IDs to HGNC gene symbols using the
# platform-specific annotation package.

gene_symbols <- mapIds(
  clariomshumantranscriptcluster.db,
  keys    = rownames(expr_matrix),
  column  = "SYMBOL",
  keytype = "PROBEID",
  multiVals = "first"       # keep one symbol per probe
)

# Attach gene symbols as a new row name
rownames(expr_matrix) <- gene_symbols

# Remove probes that could not be mapped (NA)
expr_matrix <- expr_matrix[!is.na(rownames(expr_matrix)), ]

# Remove duplicate gene symbols — keep the first occurrence
# (alternatively, collapse by mean; "first" is acceptable for exploratory work)
expr_matrix <- expr_matrix[!duplicated(rownames(expr_matrix)), ]

cat("Genes remaining after annotation:", nrow(expr_matrix), "\n")


# ============================================================
# SECTION 7 — SAMPLE METADATA & DESIGN MATRIX
# ============================================================
# 5 Vehicle controls followed by 5 MI (myocardial infarction) samples.
# Adjust the order if your CEL file names differ.

group <- factor(c(rep("Vehicle", 5), rep("MI", 5)),
                levels = c("Vehicle", "MI"))   # Vehicle is the reference
group

# Model matrix: intercept = Vehicle mean, groupMI = MI – Vehicle
design <- model.matrix(~ group)
colnames(design) <- c("Intercept", "MI_vs_Vehicle")
design


# ============================================================
# SECTION 8 — SAMPLE-LEVEL QC
# ============================================================

## 8a — Hierarchical clustering dendrogram
sample_dist <- dist(t(expr_matrix))
hc <- hclust(sample_dist, method = "complete")
plot(hc,
     main = "Sample Clustering Dendrogram",
     xlab = "",
     sub  = "",
     ylab = "Euclidean Distance")

## 8b — PCA plot (ggplot2 version)
pca     <- prcomp(t(expr_matrix), scale. = TRUE)
pca_df  <- as.data.frame(pca$x[, 1:2])
pca_df$Group  <- group
pca_df$Sample <- colnames(expr_matrix)

# Variance explained by each PC
var_explained <- round(100 * summary(pca)$importance[2, 1:2], 1)

ggplot(pca_df, aes(x = PC1, y = PC2, colour = Group, label = Sample)) +
  geom_point(size = 4) +
  geom_text_repel(size = 3) +
  scale_colour_manual(values = c("Vehicle" = "steelblue", "MI" = "firebrick")) +
  labs(
    title = "PCA — Clariom S Human (GSE233849)",
    x = paste0("PC1 (", var_explained[1], "% variance)"),
    y = paste0("PC2 (", var_explained[2], "% variance)")
  ) +
  theme_bw(base_size = 13)


# ============================================================
# SECTION 9 — DIFFERENTIAL EXPRESSION ANALYSIS (limma)
# ============================================================
# limma fits a linear model per gene and uses empirical Bayes
# moderation of variance estimates (eBayes) to improve power
# with small sample sizes.

fit  <- lmFit(expr_matrix, design)
fit2 <- eBayes(fit)

# Extract full results table (all genes, sorted by adjusted p-value)
results <- topTable(fit2,
                    coef   = "MI_vs_Vehicle",
                    number = Inf,
                    sort.by = "P")

# Attach gene symbol column for clarity
results$GeneSymbol <- rownames(results)
head(results)

# Significant genes: FDR < 5 %
sig_genes <- subset(results, adj.P.Val < 0.05)
cat("Significant DEGs (FDR < 0.05):", nrow(sig_genes), "\n")
head(sig_genes)


# ============================================================
# SECTION 10 — VISUALIZATION
# ============================================================

## 10a — MA Plot
# X-axis: average expression (A); Y-axis: log fold-change (M).
# Red points = significant genes.
plotMA(fit2, coef = "MI_vs_Vehicle",
       main = "MA Plot — MI vs Vehicle",
       status = ifelse(results$adj.P.Val < 0.05, "sig", "ns"),
       values = c("sig" = "red", "ns" = "grey70"),
       hl.pch = 16, hl.cex = 0.6)


## 10b — Volcano Plot (ggplot2, labeled top genes)
results$Significance <- "Not significant"
results$Significance[results$adj.P.Val < 0.05 & results$logFC >  1] <- "Up"
results$Significance[results$adj.P.Val < 0.05 & results$logFC < -1] <- "Down"

top_label <- head(results[order(results$adj.P.Val), ], 20)

ggplot(results, aes(x = logFC, y = -log10(adj.P.Val),
                    colour = Significance)) +
  geom_point(alpha = 0.6, size = 1.5) +
  geom_text_repel(data = top_label,
                  aes(label = GeneSymbol),
                  size = 3, max.overlaps = 20) +
  scale_colour_manual(values = c("Up"   = "firebrick",
                                  "Down" = "steelblue",
                                  "Not significant" = "grey70")) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", colour = "black") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", colour = "black") +
  labs(title = "Volcano Plot — MI vs Vehicle",
       x = "log2 Fold Change",
       y = "-log10 (Adjusted P-value)") +
  theme_bw(base_size = 13)


## 10c — Heatmap of top 50 DEGs (gene names on Y-axis)
top50 <- rownames(sig_genes)[1:min(50, nrow(sig_genes))]
top50 <- top50[top50 %in% rownames(expr_matrix)]

annotation_df <- data.frame(Group = group)
rownames(annotation_df) <- colnames(expr_matrix)

ann_colors <- list(Group = c(Vehicle = "steelblue", MI = "firebrick"))

pheatmap(
  expr_matrix[top50, ],
  scale          = "row",           # z-score per gene
  annotation_col = annotation_df,
  annotation_colors = ann_colors,
  show_rownames  = TRUE,            # gene symbols visible
  show_colnames  = TRUE,
  fontsize_row   = 8,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method        = "complete",
  main = "Top 50 DEGs — MI vs Vehicle (FDR < 0.05)"
)


# ============================================================
# SECTION 11 — SAVE RESULTS
# ============================================================

write.csv(results,   file = "DEG_results_all.csv",         row.names = TRUE)
write.csv(sig_genes, file = "DEG_results_significant.csv", row.names = TRUE)
cat("Results saved.\n")


# ============================================================
# SECTION 12 — GENE CO-EXPRESSION NETWORK
# ============================================================
# Build a simple correlation-based network from the top 50 DEGs.

## 12a — Select top genes and compute pairwise Pearson correlation
top_genes_net <- rownames(sig_genes)[1:min(50, nrow(sig_genes))]
network_data  <- expr_matrix[top_genes_net, ]
cor_matrix    <- cor(t(network_data), method = "pearson")

## 12b — Correlation heatmap
corrplot(cor_matrix,
         method  = "color",
         tl.cex  = 0.5,
         title   = "Gene Co-expression Correlation Matrix",
         mar     = c(0, 0, 2, 0))

## 12c — Threshold the correlation matrix to build adjacency
# Genes with |r| > 0.8 are considered co-expressed (connected)
adj_matrix <- abs(cor_matrix) > 0.8
diag(adj_matrix) <- FALSE   # remove self-loops

graph <- graph_from_adjacency_matrix(adj_matrix,
                                     mode     = "undirected",
                                     weighted = NULL)

# Remove isolated nodes for a cleaner plot
graph_clean <- delete.vertices(graph, degree(graph) == 0)

plot(graph_clean,
     vertex.size        = 6,
     vertex.label.cex  = 0.55,
     vertex.color       = "orange",
     vertex.frame.color = "grey40",
     edge.color         = "grey60",
     layout             = layout_with_fr,
     main               = "Gene Co-expression Network (|r| > 0.8)")

# ============================================================
# SECTION 13 — SAVE ALL PLOTS AS PNG
# ============================================================

dir.create("plots", showWarnings = FALSE)  # creates plots/ folder

# --- Raw boxplot ---
png("plots/01_raw_boxplot.png", width = 1200, height = 800, res = 150)
boxplot(data,
        main = "Raw Data — Probe Intensity Distribution",
        ylab = "log2 Intensity", las = 2, col = "lightblue")
dev.off()

# --- Normalized boxplot ---
png("plots/02_normalized_boxplot.png", width = 1200, height = 800, res = 150)
boxplot(expr_matrix,
        main = "Normalized Data — RMA",
        ylab = "log2 Intensity", las = 2, col = "lightgreen")
dev.off()

# --- Dendrogram ---
png("plots/03_dendrogram.png", width = 1200, height = 800, res = 150)
plot(hc, main = "Sample Clustering Dendrogram",
     xlab = "", sub = "", ylab = "Euclidean Distance")
dev.off()

# --- PCA ---
png("plots/04_pca.png", width = 1400, height = 1000, res = 150)
print(
  ggplot(pca_df, aes(x = PC1, y = PC2, colour = Group, label = Sample)) +
    geom_point(size = 4) +
    geom_text_repel(size = 3) +
    scale_colour_manual(values = c("Vehicle" = "steelblue", "MI" = "firebrick")) +
    labs(title = "PCA — Clariom S Human (GSE233849)",
         x = paste0("PC1 (", var_explained[1], "% variance)"),
         y = paste0("PC2 (", var_explained[2], "% variance)")) +
    theme_bw(base_size = 13)
)
dev.off()

# --- MA plot ---
png("plots/05_MA_plot.png", width = 1200, height = 800, res = 150)
plotMA(fit2, coef = "MI_vs_Vehicle",
       main = "MA Plot — MI vs Vehicle",
       status = ifelse(results$adj.P.Val < 0.05, "sig", "ns"),
       values = c("sig" = "red", "ns" = "grey70"),
       hl.pch = 16, hl.cex = 0.6)
dev.off()

# --- Volcano plot ---
png("plots/06_volcano.png", width = 1400, height = 1000, res = 150)
print(
  ggplot(results, aes(x = logFC, y = -log10(adj.P.Val), colour = Significance)) +
    geom_point(alpha = 0.6, size = 1.5) +
    geom_text_repel(data = top_label, aes(label = GeneSymbol),
                    size = 3, max.overlaps = 20) +
    scale_colour_manual(values = c("Up" = "firebrick",
                                   "Down" = "steelblue",
                                   "Not significant" = "grey70")) +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
    labs(title = "Volcano Plot — MI vs Vehicle",
         x = "log2 Fold Change", y = "-log10 (Adjusted P-value)") +
    theme_bw(base_size = 13)
)
dev.off()

# --- Heatmap ---
png("plots/07_heatmap.png", width = 1400, height = 1800, res = 150)
pheatmap(expr_matrix[top50, ],
         scale = "row",
         annotation_col = annotation_df,
         annotation_colors = ann_colors,
         show_rownames = TRUE,
         show_colnames = TRUE,
         fontsize_row = 8,
         clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean",
         clustering_method = "complete",
         main = "Top 50 DEGs — MI vs Vehicle (FDR < 0.05)")
dev.off()

# --- Correlation heatmap ---
png("plots/08_correlation_heatmap.png", width = 1400, height = 1400, res = 150)
corrplot(cor_matrix, method = "color", tl.cex = 0.5,
         title = "Gene Co-expression Correlation Matrix",
         mar = c(0, 0, 2, 0))
dev.off()

# --- Gene network ---
png("plots/09_gene_network.png", width = 2000, height = 2000, res = 200)
plot(graph_clean,
     vertex.size        = 6,
     vertex.label.cex   = 0.55,
     vertex.color       = "orange",
     vertex.frame.color = "grey40",
     edge.color         = "grey60",
     layout             = layout_with_fr,
     main               = "Gene Co-expression Network (|r| > 0.8)")
dev.off()

cat("All plots saved to plots/ folder.\n")