################################################################################
# Orlando et al. Nature 2021
# Script for Figure 1j: score CRC primary cells with signatures derived from H9 scRNAseq data and Stem Cell Index
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

# Set-up and QC

#set working directory
# downloaded "GSE81861_CRC_tumor_epithelial_cells_COUNT.csv" from GSE81861 (NCBI, GEO)
# unzipped the files using Terminal with command "gunzip file.csv.gz"
# store file in working directory

#load data
CRC_tumor <- read.csv("GSE81861_CRC_tumor_epithelial_cells_COUNT.csv", header=TRUE)

# format Gene Symbol column
Genetable <- data.frame(do.call('rbind', strsplit(as.character(CRC_tumor$X),'_',fixed=TRUE)))
CRC_tumor_2 <- CRC_tumor
CRC_tumor_2$X <- Genetable$X3
CRC_tumor_2 <- CRC_tumor_2[!duplicated(CRC_tumor_2$X), ]
Genetable <- Genetable[!duplicated(Genetable$X3), ]
# Remove decimal places from  ENSG gene names
Genetable_2 <- separate(data = Genetable, col = X3, into = c("ENSG", "right"), sep = "\\.")
# Genetable_2$X2 = gene symbols
# Genetable_2$ENSG = ENSG gene names

# put gene annotation as row name
CRC_tumor_2$X <- Genetable_2$ENSG
# remove two rows called "A"
CRC_tumor_2 <- CRC_tumor_2[!duplicated(CRC_tumor_2$X), ]
Genetable_2 <- Genetable_2[!duplicated(Genetable_2$ENSG), ]
rownames(CRC_tumor_2) <- CRC_tumor_2[,1]
CRC_tumor_2[,1] <- NULL

# succinct Genetable will become the row metadata
Genetable_2[,4:5] <- NULL
Genetable_2[,1:2] <- NULL
head(Genetable_2)
#             ENSG       X5
#1 ENSG00000000003   TSPAN6
#2 ENSG00000000005     TNMD
#3 ENSG00000000419     DPM1

# set-up stem cell experiment (sce)
all.counts <- CRC_tumor_2
all.counts <- as.matrix(all.counts)
cell_labels <- colnames(CRC_tumor_2)
genes <- Genetable_2

sce <- SingleCellExperiment(list(counts=all.counts),
                            colData=DataFrame(cell_labels),
                            rowData=DataFrame(genes),
                            metadata=list(study="GSE81861_CRC"))

# useful commands to check data:
dim(sce)
head(rownames(sce))
head(rowData(sce)$ENSG)
head(rowData(sce)$X5)
assays(sce)$counts

# QC, using mitochondrial gene list and RNA abundance
is.mito <- grepl("^MT-", rowData(sce)$X5)
summary(is.mito)
#   Mode   FALSE    TRUE 
#logical   56378      13  

sce2 <- perCellQCMetrics(sce, subsets=list(Mt=is.mito)) 
par(mfrow=c(2,2), mar=c(5.1, 4.1, 0.1, 0.1))
hist(sce2@listData[["sum"]], xlab="Library sizes", main="", 
     breaks=20, col="grey80", ylab="Number of cells")
hist(sce2@listData[["detected"]], xlab="Number of expressed genes", main="", 
     breaks=20, col="grey80", ylab="Number of cells")
hist(sce2@listData[["subsets_Mt_percent"]], xlab="Mitochondrial proportion (%)", 
     ylab="Number of cells", breaks=20, main="", col="grey80")
# Saved as CRC_01_QCmetrics.pdf
ggsave("./CRC_01_QCmetrics.pdf", width = 5, height = 4)

libsize.drop <- isOutlier(sce2@listData[["sum"]], nmads=5, type="lower", log=TRUE)
feature.drop <- isOutlier(sce2@listData[["detected"]], nmads=5, type="lower", log=TRUE)
mito.drop <- isOutlier(sce2@listData[["subsets_Mt_percent"]], nmads=5, type="higher", log=TRUE)
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
#0.01094 0.28421 0.70205 1.00000 1.52100 4.84494
plot(sizeFactors(sce3), sce$total_counts, log="xy",
     ylab="Library size (thousands)", xlab="Size factor")
sce3 <- logNormCounts(sce3)

# Modeling and removing technical noise 
var.out <- modelGeneVar(sce3)
plot(var.out$mean, var.out$total, xlab="Mean log-expression", ylab="Variance of log-expression")
curve(metadata(var.out)$trend(x), col="blue", add=TRUE)
# Saved as CRC_02_variance_of_normalized_log-expression_values_forEachGene.pdf
ggsave("./CRC_02_variance_of_normalized_log-expression_values_forEachGene.pdf", width = 5, height = 4)

# Get the top 10% of genes.
top.hvgs <- getTopHVGs(var.out, prop=0.1)
# Get all genes with positive biological components.
top.hvgs3 <- getTopHVGs(var.out, var.threshold=0)

# Automated PC choice
dec2 <- var.out
sced <- denoisePCA(sce3, dec2, subset.row=getTopHVGs(dec2, prop=0.1))
ncol(reducedDim(sced, "PCA"))
# [1] 5

# Assigning the colLabels
rowData(sce)
#DataFrame with 56389 rows and 2 columns
#                           ENSG          X5
#                    <character> <character>
#ENSG00000000003 ENSG00000000003      TSPAN6
rowData(sced)$ENSEMBL <- rownames(sced)
rownames(sced) <- rowData(sced)$X5
head(rownames(sced))
# [1] "TSPAN6"   "TNMD"     "DPM1"     "SCYL3"    "C1orf112" "FGR"  

g <- buildSNNGraph(sced, use.dimred="PCA")
cluster <- igraph::cluster_walktrap(g)$membership
colLabels(sced) <- factor(cluster)
table(colLabels(sced))
#  1  2  3  4  5  6 
# 53 93 57 26 24 19
sced$cluster_igraph <- factor(cluster)
sced <- runTSNE(sced, dimred="PCA")

# plot PCA
plotPCA(sced, colour_by="cluster_igraph", text_by="cluster_igraph")

# plot t-SNE:
temp <- as.data.frame(reducedDim(sced, "TSNE"))
ggplot(temp, aes(x=V1, y=V2)) + geom_point(aes(colour=sced$cluster_igraph), size=1)+theme_test()

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
write.table(temp, file = "CRC_logcounts.txt", sep="\t", row.names=TRUE, col.names=TRUE)
# Manually added ["gene_id" and (tab)] at the head of the "CRC_logcounts.txt" file 

# Manually added a header with gene_id, switched header=true for the "pcbc-stemsig.tsv"
w <- read.delim("pcbc-stemsig.tsv", header=TRUE, row.names=1) %>% as.matrix() %>% drop()
f <- function( v ) unlist( lapply( strsplit( v, "\\|" ), "[[", 1 ) )
X <- read.delim("CRC_logcounts.txt", as.is=TRUE, check.names=FALSE, sep="\t") %>% drop() %>% 
  filter( !grepl( "\\?", gene_id ) ) %>%      ## Drop genes with no mapping to HUGO
  mutate( gene_id = f( gene_id ) ) %>%        ## Clip gene ids to HUGO
  filter( gene_id %in% names(w) )         ## Reduce to the signature's gene set

rownames(X) <- NULL
X <- X %>% tibble::column_to_rownames( "gene_id" ) %>% as.matrix()
stopifnot( all( rownames(X) %in% names(w) ) )
w <- w[ rownames(X) ]

# Score via Spearman correlation, and scale between 0 and 1
s <- apply( X, 2, function(z) {cor( z, w, method = "sp", use = "complete.obs" )} )
s <- s - min(s)
s <- s / max(s)

write.table(cbind(s), file = "CRC_stemscore.tsv", sep = "\t", quote = FALSE, col.names = FALSE)

mRNAsi <- read.delim("CRC_stemscore.tsv", header=FALSE, row.names=1)
cellID <- colnames(sced)
mRNAsi <- mRNAsi$V2[match(cellID, rownames(mRNAsi))]
mRNAsi <- t(mRNAsi)
colData(sced)$mRNAsi <- mRNAsi

# Plot stem cell index (mRNAsi) using t-SNE reduction
temp <- as.data.frame(reducedDim(sced, "TSNE"))
ggplot(temp, aes(x=V1, y=V2)) + geom_point(aes(colour=sced$mRNAsi), size=1)+theme_test() +
    scale_colour_gradient2(low = "blue", mid = "lightgrey", high = "red", 
    midpoint = 0.5, space = "Lab",guide = "colourbar", aesthetics = "colour")
ggsave("./CRC_03_tSNE_mRNAsi.pdf", width =5, height = 4)
# This plot is shown in Orlando et al. in Figure 1j, on the left

################################################################################

# Original analysis: score CRC primary cells with signatures derived from H9 scRNAseq data

# Uncommitted score, based on clusters 1-4
read.csv("List_uncommitted.csv", header=TRUE)
Primitive1 = readLines("List_uncommitted.csv")
Primitive1.genes <- list(Primitive1,rownames(sced))
Primitive1.genes <- Reduce(f = intersect, x = Primitive1.genes)
Primitive1.list <- as.list(Primitive1.genes)
sced_Seurat <- AddModuleScore(object=sced_Seurat, list(features=Primitive1.list), ctrl = 5, name = 'Primitive1_')
FeatureScatter(sced_Seurat, 'Primitive1_1', 'mRNAsi', cols='orange')

# Add regression line 
temp2 <- FetchData(sced_Seurat, c("Primitive1_1", "mRNAsi"))
temp2 <- as.data.frame(temp2)
x <- temp2$Primitive1_1
y <- temp2$mRNAsi
plot(x, y, main = "CRC",xlab = "Primitive1_1", ylab = "mRNAsi",
     pch = 19, frame = FALSE, col = "red", xlim=c(0,0.25))
abline(lm(y ~ x, data = temp2), col = "black")

# Pearson's product-moment correlation
cor.test(sced_Seurat$Primitive_1, sced_Seurat$mRNAsi, method=c("pearson"))

# Saved as "CRC_04_correlation_mRNAsi_uncommitted_regression.pdf"
# This plot is shown in Orlando et al. in Figure 1j, in the centre

# Primed score, based on clusters 6-10
read.csv("List_primed.csv", header=TRUE)
Diff = readLines("List_primed.csv")
Diff.genes <- list(Diff,rownames(sced))
Diff.genes <- Reduce(f = intersect, x = Diff.genes)
Diff.list <- as.list(Diff.genes)
sced_Seurat <- AddModuleScore(object=sced_Seurat, list(features=Diff.list), ctrl = 5, name = 'Diff_')
FeatureScatter(sced_Seurat, 'Diff_1', 'mRNAsi', cols='orange')

# Add regression line
temp2 <- FetchData(sced_Seurat, c("Diff_1", "mRNAsi"))
temp2 <- as.data.frame(temp2)
x <- temp2$Diff_1
y <- temp2$mRNAsi
plot(x, y, main = "CRC",xlab = "Diff_1", ylab = "mRNAsi",
     pch = 19, frame = FALSE, col = "red", xlim=c(-0.06,0.1))
abline(lm(y ~ x, data = temp2), col = "black")

# Pearson's product-moment correlation
cor.test(sced_Seurat$Diff_1, sced_Seurat$mRNAsi, method=c("pearson"))

# Saved as "CRC_05_correlation_mRNAsi_primed_regression.pdf"
# This plot is shown in Orlando et al. in Figure 1j, on the right

# End of Figure 1j