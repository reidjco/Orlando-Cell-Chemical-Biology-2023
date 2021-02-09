################################################################################
# Orlando et al. Nature 2021
# Script for Figure 1i: score GBM primary cells with signatures derived from H9 scRNAseq data and Stem Cell Index
# written by Dr. Jennifer Reid
################################################################################

# only install the following packages once:
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("DropletUtils")
BiocManager::install("SingleCellExperiment")
BiocManager::install("scater")
BiocManager::install("scran")
BiocManager::install("AnnotationDbi")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("pheatmap")
BiocManager::install("gelnet")
BiocManager::install("dplyr")
BiocManager::install("gdata")
BiocManager::install("ggplot")
BiocManager::install("ggpubr")
BiocManager::install("Seurat")

# load libraries:
library(SingleCellExperiment)
library(DropletUtils)
library(scater)
library(scran)
library(tidyr)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(pheatmap)
library(gelnet)
library(dplyr)
library(gdata)
library(ggplot2)
library(Seurat)
library(ggpubr)

################################################################################

# Sample Extraction

# set working directory
# download GSE141946_scRNAseq_counts.txt from GSE141946 (NCBI, GEO)
# unzipped the files using Terminal with "gunzip GSE141946_scRNAseq_counts.txt.gz" and "gunzip GSE141946_scRNAseq_metadata.csv.gz"
# store file in working directory

# load data
GBM_tumor <- read.table("GSE141946_scRNAseq_counts.txt", header=TRUE,sep = '\t')
# this file is very large, so one primary sample was selected

# subset GBM tumor # 8165C
GBM_8165C <- GBM_tumor[,9456:12684]
write.csv(GBM_8165C, file = "GBM_8165C_subset.csv")
GBM2 <- GBM_8165C

GBM2 <- read.csv("GBM_8165C_subset.csv", header=TRUE)

# Check extraction:
# top looks good:
GBM2[0:3,0:3]
#         X8165C.Tumor.1 X8165C.Tumor.2 X8165C.Tumor.3
#A1BG                  0              0              0
#A1BG-AS1              0              0              0
#A2M  

# bottom looks good:
dim(GBM2)
#[1] 40604  3229
GBM2[40602:40604, 3226:3229]
#       X8165C.Tumor.3225 X8165C.Tumor.3226 X8165C.Tumor.3227 X8165C.Tumor.3228
# 40602                 0                 0                 0                 0
# 40603                 0                 0                 0                 0
# 40604                 0                 0                 0                 0

################################################################################

# Set-up and QC

# gene symbol is the first column
rownames(GBM2) <- GBM2[,1]
GBM2[,1] <- NULL

# remove any potentially duplicated rows
GBM2 <- GBM2[!duplicated(rownames(GBM2)), ]

# set-up stem cell experiment (sce)
all.counts <- GBM2
all.counts <- as.matrix(all.counts)
cell_labels <- colnames(GBM2)
genes <- rownames(GBM2)

sce <- SingleCellExperiment(list(counts=all.counts),
                            colData=DataFrame(cell_labels),
                            rowData=DataFrame(genes),
                            metadata=list(study="GSE131928_GBM2"))

# useful commands to check data:
dim(sce)
colData(sce)
assays(sce)$counts

# to get gene symbols:
head(rowData(sce))

# QC, using mitochondrial gene list and RNA abundance
is.mito <- grepl("^MT-", rownames(sce))
summary(is.mito)
#   Mode   FALSE    TRUE 
#logical   40573      31 

sce2 <- perCellQCMetrics(sce, subsets=list(Mt=is.mito)) 
par(mfrow=c(2,2), mar=c(5.1, 4.1, 0.1, 0.1))
hist(sce2@listData[["sum"]], xlab="Library sizes", main="", 
     breaks=20, col="grey80", ylab="Number of cells")
hist(sce2@listData[["detected"]], xlab="Number of expressed genes", main="", 
     breaks=20, col="grey80", ylab="Number of cells")
hist(sce2@listData[["subsets_Mt_percent"]], xlab="Mitochondrial proportion (%)", 
     ylab="Number of cells", breaks=20, main="", col="grey80")
# Saved as GBM_01_QCmetrics.pdf
ggsave("GBM_01_QCmetrics.pdf", width = 5, height = 4)

libsize.drop <- isOutlier(sce2@listData[["sum"]], nmads=3, type="lower", log=TRUE)
feature.drop <- isOutlier(sce2@listData[["detected"]], nmads=3, type="lower", log=TRUE)
mito.drop <- isOutlier(sce2@listData[["subsets_Mt_percent"]], nmads=3, type="higher", log=TRUE)
summary(mito.drop)
summary(feature.drop)
summary(libsize.drop)
# all are FALSE and is in agreement the QC performed by the original authors of this study

# Normalization 
sce3 <- sce
clusters <- quickCluster(sce3, min.mean=0.1, method="igraph")
sce3 <- computeSumFactors(sce3, cluster=clusters, min.mean=0.1)
summary(sizeFactors(sce3))
#    Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.2045  0.3897  0.6058  1.0000  1.1235  7.3582 
plot(sizeFactors(sce3), sce$total_counts, log="xy",
     ylab="Library size (thousands)", xlab="Size factor")
sce3 <- logNormCounts(sce3)

# Modeling and removing technical noise 
var.out <- modelGeneVar(sce3)
plot(var.out$mean, var.out$total, xlab="Mean log-expression", ylab="Variance of log-expression")
curve(metadata(var.out)$trend(x), col="blue", add=TRUE)
# Saved as GBM_02_variance_of_normalized_log-expression_values_forEachGene.pdf
ggsave("GBM_02_variance_of_normalized_log-expression_values_forEachGene.pdf", width = 5, height = 4)

# Get the top 10% of genes.
top.hvgs <- getTopHVGs(var.out, prop=0.1)
# Get all genes with positive biological components.
top.hvgs3 <- getTopHVGs(var.out, var.threshold=0)

# Automated PC choice
dec2 <- var.out
sced <- denoisePCA(sce3, dec2, subset.row=getTopHVGs(dec2, prop=0.1))
ncol(reducedDim(sced, "PCA"))
# [1] 6

# Assigning the colLabels
g <- buildSNNGraph(sced, use.dimred="PCA", k=50)
cluster <- igraph::cluster_walktrap(g)$membership
colLabels(sced) <- factor(cluster)
table(colLabels(sced))
# 1   2   3   4   5   6 
# 571 736 881 103 664 274
sced$cluster_igraph <- factor(cluster)
sced <- runTSNE(sced, dimred="PCA")

# plot PCA
plotPCA(sced, colour_by="cluster_igraph", text_by="cluster_igraph")

# plot t-SNE:
temp <- as.data.frame(reducedDim(sced, "TSNE"))
ggplot(temp, aes(x=V1, y=V2)) + geom_point(aes(colour=sced$cluster_igraph), size=0.5)+theme_test()

# Draw modularity-heatmaps for the 2 methods.
ratio <- clusterModularity(g, cluster, as.ratio=TRUE)
pheatmap(log10(ratio+1), cluster_cols=FALSE, cluster_rows=FALSE,
         col=rev(heat.colors(100)))

################################################################################

# Calculate Stem Cell Index
# This section is adapted from Malta et al. Cell. 2018 Apr 5;173(2):338-354.e15. doi: 10.1016/j.cell.2018.03.034. 

temp <- logcounts(sced)
temp <- as.matrix(temp)
temp <- temp[!duplicated(row.names(temp)), ]
write.table(temp, file = "GBM_logcounts.txt", sep="\t", row.names=TRUE, col.names=TRUE)
# Manually added ["gene_id" and (tab)] at the head of the "GBM_logcounts.txt" file 

# Manually added a header with gene_id, switched header=true for the "pcbc-stemsig.tsv"
w <- read.delim("pcbc-stemsig.tsv", header=TRUE, row.names=1) %>% as.matrix() %>% drop()
f <- function( v ) unlist( lapply( strsplit( v, "\\|" ), "[[", 1 ) )
X <- read.delim("GBM_logcounts.txt", as.is=TRUE, check.names=FALSE, sep="\t") %>% drop() %>% 
  filter( !grepl( "\\?", gene_id ) ) %>%      
  mutate( gene_id = f( gene_id ) ) %>%        
  filter( gene_id %in% names(w) )         

rownames(X) <- NULL
X <- X %>% tibble::column_to_rownames( "gene_id" ) %>% as.matrix()
stopifnot( all( rownames(X) %in% names(w) ) )
w <- w[ rownames(X) ]

# Score via Spearman correlation, and scale between 0 and 1
s <- apply( X, 2, function(z) {cor( z, w, method = "sp", use = "complete.obs" )} )
s <- s - min(s)
s <- s / max(s)

write.table(cbind(s), file = "GBM_stemscore.tsv", sep = "\t", quote = FALSE, col.names = FALSE)

mRNAsi <- read.delim("GBM_stemscore.tsv", header=FALSE, row.names=1)
cellID <- colnames(sced)
mRNAsi <- mRNAsi$V2[match(cellID, rownames(mRNAsi))]
mRNAsi <- t(mRNAsi)
colData(sced)$mRNAsi <- mRNAsi

# Plot stem cell index (mRNAsi) using t-SNE reduction
temp <- as.data.frame(reducedDim(sced, "TSNE"))
ggplot(temp, aes(x=V1, y=V2)) + geom_point(aes(colour=sced$mRNAsi), size=0.5)+theme_test() +
  scale_colour_gradient2(low = "blue", mid = "lightgrey", high = "red", 
                         midpoint = 0.5, space = "Lab",guide = "colourbar", aesthetics = "colour")
ggsave("GBM_03_tSNE_mRNAsi.pdf", width =5, height = 4)
# This plot is shown in Orlando et al. in Figure 1i, on the left

################################################################################

# Original analysis: score AML primary cells with signatures derived from H9 scRNAseq data

# Uncommitted score, based on clusters 1-4
read.csv("List_uncommitted.csv", header=TRUE)
Primitive1 = readLines("List_uncommitted.csv")
Primitive1.genes <- list(Primitive1,rownames(sced))
Primitive1.genes <- Reduce(f = intersect, x = Primitive1.genes)
Primitive1.list <- as.list(Primitive1.genes)
sced_Seurat <- AddModuleScore(object=sced_Seurat, list(features=Primitive1.list), ctrl = 5, name = 'Primitive1_')
FeatureScatter(sced_Seurat, 'Primitive1_1', 'mRNAsi', cols='green4')

# Add regression line
temp2 <- FetchData(sced_Seurat, c("Primitive1_1", "mRNAsi"))
temp2 <- as.data.frame(temp2)
x <- temp2$Primitive1_1
y <- temp2$mRNAsi
plot(x, y, main = "GBM",xlab = "Primitive1_1", ylab = "mRNAsi",
     pch = 19, frame = FALSE, col = "green4", xlim=c(-0.05,0.15))
abline(lm(y ~ x, data = temp2), col = "black")

# Pearson's product-moment correlation
cor.test(sced_Seurat$Primitive1_1, sced_Seurat$mRNAsi, method=c("pearson"))

# Saved as "GBM_04_correlation_mRNAsi_uncommitted_regression.pdf"
# This plot is shown in Orlando et al. in Figure 1i, in the centre

# Primed score, based on clusters 6-10
read.csv("List_primed.csv", header=TRUE)
Diff = readLines("List_primed.csv")
Diff.genes <- list(Diff,rownames(sced))
Diff.genes <- Reduce(f = intersect, x = Diff.genes)
Diff.list <- as.list(Diff.genes)
sced_Seurat <- AddModuleScore(object=sced_Seurat, list(features=Diff.list), ctrl = 5, name = 'Diff_')
FeatureScatter(sced_Seurat, 'Diff_1', 'mRNAsi', cols='green4')

# Add regression line to the above plot
temp2 <- FetchData(sced_Seurat, c("Diff_1", "mRNAsi"))
temp2 <- as.data.frame(temp2)
x <- temp2$Diff_1
y <- temp2$mRNAsi
plot(x, y, main = "GBM",xlab = "Diff_1", ylab = "mRNAsi",
     pch = 19, frame = FALSE, col = "green4")
abline(lm(y ~ x, data = temp2), col = "black")

# Pearson's product-moment correlation
cor.test(sced_Seurat$Diff_1, sced_Seurat$mRNAsi, method=c("pearson"))

# Saved as "GBM_05_correlation_mRNAsi_primed_regression.pdf"
# This plot is shown in Orlando et al. in Figure 1i, on the right

# End of Figure 1i