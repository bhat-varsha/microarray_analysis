#==============================
# INSTALL REQUIRED PACKAGES
# ==============================
install.packages("BiocManager")
BiocManager::install("affy")
BiocManager::install("gcrma")
BiocManager::install("limma")
BiocManager::install("oligo")
BiocManager::install("pd.clariom.s.human")
BiocManager::install("clariomshumantranscriptcluster.db")
install.packages("pheatmap")
# ==============================
# LOAD LIBRARIES
# ==============================
library(affy)
library(gcrma)
library(limma)
library(oligo)
library(pd.clariom.s.human)
library(clariomshumantranscriptcluster.db)
library(pheatmap)
# ==============================
# SET WORKING DIRECTORY
# ==============================
getwd()
setwd("GSE233849_RAW")
list.files()
# ==============================
# READ CEL FILES
# ==============================
data <- read.celfiles(list.celfiles())
data
# Boxplot before normalization
boxplot(data,main="Raw Data Boxplot",las=2)
# ==============================
# NORMALIZATION (RMA)
# ==============================
norm_data <- rma(data)
norm_data
# ==============================
# EXTRACT EXPRESSION MATRIX
# ==============================
expr_matrix <- exprs(norm_data)
dim(expr_matrix)
head(expr_matrix)
# ==============================
# BOXPLOT AFTER NORMALIZATION
# ==============================
boxplot(expr_matrix,main="Normalized Data Boxplot",las=2)
# ==============================
# SAMPLE CLUSTERING
# ==============================
# Calculate sample distances
sample_dist <- dist(t(expr_matrix))
# Hierarchical clustering
hc <- hclust(sample_dist)
# Plot dendrogram
plot(hc,main="Sample Clustering Dendrogram",xlab="Samples")
# ==============================
# PCA ANALYSIS
# ==============================
pca <- prcomp(t(expr_matrix),
              scale.=TRUE)
plot(pca$x[,1],
     pca$x[,2],
     col=c(rep("blue",5), rep("red",5)),
     pch=19,
     xlab="PC1",
     ylab="PC2",
     main="PCA Plot")
legend("topright",
       legend=c("Vehicle","MI"),
       col=c("blue","red"),
       pch=19)
# ==============================
# CREATE GROUPS
# ==============================
group <- c("Vehicle","Vehicle","Vehicle","Vehicle","Vehicle",
           "MI","MI","MI","MI","MI")
group <- factor(group)
group
# ==============================
# DESIGN MATRIX
# ==============================
design <- model.matrix(~group)
design
# ==============================
# LIMMA DIFFERENTIAL EXPRESSION
# ==============================
fit <- lmFit(expr_matrix, design)
fit2 <- eBayes(fit)
# ==============================
# GET DIFFERENTIAL EXPRESSION RESULTS
# ==============================
results <- topTable(fit2,
                    coef=2,
                    number=Inf)
head(results)
# ==============================
# SIGNIFICANT GENES
# ==============================
sig_genes <- results[results$adj.P.Val < 0.05, ]
nrow(sig_genes)
head(sig_genes)
# ==============================
# MA PLOT
# ==============================
plotMA(fit2)
# ==============================
# VOLCANO PLOT
# ==============================
volcanoplot(fit2,
            coef=2,
            main="Volcano Plot")
# ==============================
# HEATMAP OF TOP GENES
# ==============================
top_genes <- rownames(sig_genes)[1:min(50, nrow(sig_genes))]
top_genes <- top_genes[top_genes %in% rownames(expr_matrix)]
length(top_genes)
top_genes
nrow(sig_genes)
annotation_df <- data.frame(Group=group)
rownames(annotation_df) <- colnames(expr_matrix)
pheatmap(expr_matrix[top_genes, ],
         scale="row",
         annotation_col=annotation_df)
# ==============================
# ANNOTATE GENE SYMBOLS
# ==============================
library(AnnotationDbi)
gene_symbols <- mapIds(
  clariomshumantranscriptcluster.db,
  keys=rownames(results),
  column="SYMBOL",
  keytype="PROBEID",
  multiVals="first"
)
results$GeneSymbol <- gene_symbols
head(results)
# ==============================
# SAVE RESULTS
# ==============================
write.csv(results,
          "DEG_results.csv")
write.csv(sig_genes,
          "Significant_Genes.csv")
#___________________________________________________
#Step 1 — Install Packages
install.packages("igraph")
install.packages("corrplot")
BiocManager::install("WGCNA")
#Step 2 — Load Libraries
library(WGCNA)
library(igraph)
library(corrplot)
#Step 3 — Select Top Significant Genes
#Use top 50 or 100 genes.
top_genes <- rownames(sig_genes)[1:50]
network_data <- expr_matrix[top_genes, ]
#Step 4 — Calculate Correlation Matrix
cor_matrix <- cor(t(network_data))
head(cor_matrix)
#This checks which genes behave similarly.
#Step 5 — Visualize Correlation Heatmap
corrplot(cor_matrix,method="color",tl.cex=0.5)
#This is your first simple gene network visualization.
#Step 6 — Build Gene Network
adjacency_matrix <- cor_matrix > 0.8
graph <- graph_from_adjacency_matrix(adjacency_matrix,
                                     mode="undirected")
plot(graph,
     vertex.size=5,
     vertex.label.cex=0.6)


top_genes <- rownames(sig_genes)[1:50]

top_genes <- top_genes[!is.na(top_genes)]

top_genes <- top_genes[top_genes %in% rownames(expr_matrix)]

length(top_genes)

heatmap_data <- expr_matrix[top_genes, ]

dim(heatmap_data)

sum(is.na(heatmap_data))

heatmap_data <- heatmap_data[complete.cases(heatmap_data), ]

dim(heatmap_data)

annotation_df <- data.frame(Group=group)
rownames(annotation_df) <- colnames(heatmap_data)

pheatmap(
  heatmap_data,
  scale="row",
  annotation_col=annotation_df,
  main="Heatmap of Differentially Expressed Genes"
)

# ==============================
# GENE ANNOTATION
# ==============================
library(AnnotationDbi)
library(clariomshumantranscriptcluster.db)
# Map probe IDs to gene symbols
gene_symbols <- mapIds(
  clariomshumantranscriptcluster.db,
  keys=rownames(expr_matrix),
  column="SYMBOL",
  keytype="PROBEID",
  multiVals="first"
)

# Replace probe IDs with gene symbols
rownames(expr_matrix) <- gene_symbols

# Remove missing gene names
expr_matrix <- expr_matrix[!is.na(rownames(expr_matrix)), ]

# Remove duplicate genes
expr_matrix <- expr_matrix[!duplicated(rownames(expr_matrix)), ]

 # ==============================
# RERUN LIMMA
# ==============================

fit <- lmFit(expr_matrix, design)

fit2 <- eBayes(fit)

results <- topTable(fit2,
                    coef=2,
                    number=Inf)

sig_genes <- results[results$adj.P.Val < 0.05, ]

# ==============================
# HEATMAP
# ==============================

top_genes <- rownames(sig_genes)[1:min(50, nrow(sig_genes))]

top_genes <- top_genes[top_genes %in% rownames(expr_matrix)]

annotation_df <- data.frame(Group=group)

rownames(annotation_df) <- colnames(expr_matrix)

pheatmap(expr_matrix[top_genes, ],
         scale="row",
         annotation_col=annotation_df,
         main="Heatmap of Differentially Expressed Genes")





#other parts

install.packages("BiocManager")
BiocManager::install("affy")
BiocManager::install("gcrma")
BiocManager::install("limma")
library(affy)
library(gcrma)
library(limma)
setwd("D:/msc/bioinfo/genomics/practice/GSE233849_RAW")
list.files()
data <- ReadAffy()
data
BiocManager::install("clariomshumantranscriptcluster.db")
BiocManager::install("pd.clariom.s.human")
library(pd.clariom.s.human)
BiocManager::install("oligo")
library(oligo)
data <- read.celfiles(list.celfiles())
data

norm_data <- rma(data)
norm_data

expr_matrix <- exprs(norm_data)
expr_matrix
dim(expr_matrix)
head(expr_matrix)
group <- c("Vehicle","Vehicle","Vehicle","Vehicle","Vehicle",
           "MI","MI","MI","MI","MI")
group
group <- factor(group)
group
design <- model.matrix(~group)
design
fit <- lmFit(expr_matrix, design)
fit2 <- eBayes(fit)
results <- topTable(fit2, coef=2, number=Inf)
head(results)
sig_genes <- results[results$adj.P.Val < 0.05, ]
nrow(sig_genes)
volcanoplot(fit2, coef=2)
write.csv(results, "DEG_results.csv")
head(sig_genes)
install.packages("pheatmap")
library(pheatmap)
top_genes <- rownames(sig_genes)[1:50]
pheatmap(expr_matrix[top_genes, ])

#ml model training
install.packages("randomForest")
library(randomForest)
top50 <- rownames(sig_genes)[1:50]
ml_data <- t(expr_matrix[top50, ])
ml_data <- as.data.frame(ml_data)
ml_data$group <- group
head(ml_data)
colnames(ml_data) <- make.names(colnames(ml_data))
rf_model <- randomForest(group ~ ., data=ml_data)
rf_model
predictions <- predict(rf_model)
table(predictions, group)
mean(predictions == group)
importance(rf_model)
varImpPlot(rf_model)
