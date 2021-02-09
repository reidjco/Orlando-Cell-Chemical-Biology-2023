################################################################################
# Orlando et al. Nature 2021
# R script used for merged primary human AML and H9 scRNAseq data analysis 
# Data panel shown in Extended Data Figure 3f
# written by: Dr. Mio Nakanishi
# Date: October 3, 2019
################################################################################

# The AML count matrix file is publicly available on NCBI, GEO accession: GSE132567

# The H9 count matrix file is publicly available on NCBI, GEO accession: GSE126022 
# Originally published: Nakanishi et al. Cell 2019. PMID: 30982595

################################################################################

# load libraries
library(DropletUtils)
library(scran)
library(scater)

################################################################################

## LOADING THE COUNT MATRIX AND SEPARATION OF BULK H9 CELLS.
## These steps were done according to the script used for Figure 1 H9 analysis

## After the loading/QC:
assignments.H9 <- assignments
cellcodes.H9 <- cellcodes
feature.drop.H9 <- feature.drop
hs.pairs.H9 <- hs.pairs
is.mito.H9 <- is.mito
libsize.drop.H9 <- libsize.drop
mito.drop.H9 <- mito.drop
sce.H9 <- sce

rm(assignments)
rm(cellcodes)
rm(feature.drop)
rm(hs.pairs)
rm(is.mito)
rm(libsize.drop)
rm(mito.drop)
rm(sce)

ave.counts.H9 <- calcAverage(sce.H9, use_size_factors=FALSE)

rowData(sce.H9)$ave.counts.H9 <- ave.counts.H9
to.keep.H9 <- ave.counts.H9 > 0
sce.H9 <- sce.H9[to.keep.H9,]
summary(to.keep.H9)

## computeSumFactors
clusters.H9 <- quickCluster(sce.H9, min.mean=0.1, method="igraph")
sce.H9 <- computeSumFactors(sce.H9, cluster=clusters.H9, min.mean=0.1)
summary(sizeFactors(sce.H9))

sce.H9 <- normalize(sce.H9)

## Correcting batch effects
var.fit.nospike.H9 <- trendVar(sce.H9, parametric=TRUE, use.spikes=FALSE, loess.args=list(span=0.2))
var.dec.nospike.H9 <- decomposeVar(sce.H9, var.fit.nospike.H9)

var.dec.nospike.H9$Symbol <- rowData(sce.H9)$Symbol
var.dec.nospike.H9 <- var.dec.nospike.H9[order(var.dec.nospike.H9$bio, decreasing=TRUE),]
head(var.dec.nospike.H9)

################################################################################

## Load count matrix for the primary human AML sample, called 'A241_Dx'
sce.A241Dx <- read10xCounts("~/A241Dx_counts.txt", col.names = TRUE)
dim(sce.A241Dx)

## Making mitochondrial gene list
is.mito.A241Dx <- grepl("^hg19_MT-", rowData(sce.A241Dx)$Symbol)
summary(is.mito.A241Dx)

sce.A241Dx <- calculateQCMetrics(sce.A241Dx, feature_controls=list(Mt=is.mito.A241Dx))

libsize.drop.A241Dx <- isOutlier(sce.A241Dx$total_counts, nmads=3, type="lower", log=TRUE)
feature.drop.A241Dx <- isOutlier(sce.A241Dx$total_features_by_counts, nmads=3, type="lower", log=TRUE)
mito.drop.A241Dx <- isOutlier(sce.A241Dx$pct_counts_Mt, nmads=3, type="higher", log=TRUE)
sce.A241Dx <- sce.A241Dx[,!(libsize.drop.A241Dx | feature.drop.A241Dx | mito.drop.A241Dx)]
data.frame(ByLibSize=sum(libsize.drop.A241Dx), ByFeature=sum(feature.drop.A241Dx), ByMito=sum(mito.drop.A241Dx), Remaining=ncol(sce.A241Dx))

ensembl.A241Dx <- as.data.frame(rowData(sce.A241Dx)$ID)
colnames(ensembl.A241Dx) <- "ENSEMBL"
rownames(ensembl.A241Dx) <- ensembl.A241Dx$ENSEMBL
ensembl.A241Dx$ENSEMBL <- as.factor(gsub(pattern="hg19_", replacement="", ensembl.A241Dx$ENSEMBL))
head(ensembl.A241Dx$ENSEMBL)
rowData(sce.A241Dx)$ENSEMBL <- ensembl.A241Dx$ENSEMBL

hs.pairs.A241Dx <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
assignments.A241Dx <- cyclone(sce.A241Dx, hs.pairs.A241Dx, gene.names=rowData(sce.A241Dx)$ENSEMBL)
sce.A241Dx$phases <- assignments.A241Dx$phases
table(sce.A241Dx$phases)

ave.counts.A241Dx <- calcAverage(sce.A241Dx, use_size_factors=FALSE)
rowData(sce.A241Dx)$ave.counts.A241Dx <- ave.counts.A241Dx
to.keep.A241Dx <- ave.counts.A241Dx > 0
sce.A241Dx <- sce.A241Dx[to.keep.A241Dx,]
summary(to.keep.A241Dx)

clusters.A241Dx <- quickCluster(sce.A241Dx, min.mean=0.1, method="igraph")
sce.A241Dx <- computeSumFactors(sce.A241Dx, cluster=clusters.A241Dx, min.mean=0.1)
summary(sizeFactors(sce.A241Dx))

sce.A241Dx <- normalize(sce.A241Dx)

var.fit.nospike.A241Dx <- trendVar(sce.A241Dx, parametric=TRUE, use.spikes=FALSE, loess.args=list(span=0.2))
var.dec.nospike.A241Dx <- decomposeVar(sce.A241Dx, var.fit.nospike.A241Dx)

var.dec.nospike.A241Dx$Symbol <- rowData(sce.A241Dx)$Symbol
var.dec.nospike.A241Dx <- var.dec.nospike.A241Dx[order(var.dec.nospike.A241Dx$bio, decreasing=TRUE),]
head(var.dec.nospike.A241Dx)

## If necessary:
rownames(var.dec.nospike.A241Dx) <- gsub(pattern="hg19_", replacement="", rownames(var.dec.nospike.A241Dx))
rownames(sce.A241Dx) <- gsub(pattern="hg19_", replacement="", rownames(sce.A241Dx))

################################################################################

## FEATURE SELECTION ACROSS BATCHES

universe <- intersect(rownames(var.dec.nospike.H9), rownames(var.dec.nospike.A241Dx))
mean.bio <- (var.dec.nospike.H9[universe,"bio"] + var.dec.nospike.A241Dx[universe,"bio"])/2
chosen <- universe[mean.bio > 0]
length(chosen)

rescaled <- multiBatchNorm(sce.H9[universe,], sce.A241Dx[universe,])
rescaled.H9 <- rescaled[[1]]
rescaled.A241Dx <- rescaled[[2]]

## Perform MNN-based correction
set.seed(100) 
original <- list(
  H9=logcounts(rescaled.H9)[chosen,],
  A241Dx=logcounts(rescaled.A241Dx)[chosen,]
)

mnn.out <- do.call(fastMNN, c(original, list(k=20, d=50, approximate=TRUE)))
dim(mnn.out$corrected)

omat <- do.call(cbind, original)
sce <- SingleCellExperiment(list(logcounts=omat))
reducedDim(sce, "MNN") <- mnn.out$corrected
sce$Batch <- as.character(mnn.out$batch)

set.seed(100)
csce <- runTSNE(sce, use_dimred="MNN", perplexity=50)

A241.list <- read.table(file="A241_LSC17_information.tsv", header=FALSE, sep="\t")

temp <- as.data.frame(reducedDim(csce, "TSNE"))
temp$Batch <- colData(csce)$Batch
temp$pos <- A241.list$V3[match(colnames(csce), A241.list$V1)]
temp$pos <- as.character(temp$pos)
temp$pos[temp$Batch == "H9"] <- "ESC"

png("BulkH9_A241Dx_tSNE_LSC17pos.png", width = 1000, height = 1000)
ggplot(temp, aes(x=V1, y=V2)) + geom_point(size=3, aes(color=pos)) + theme_classic(base_size = 36) + scale_color_manual(limits=c("TRUE", "FALSE","ESC"), values=c("#F2519B", "#B60E5B", "snow2"), labels=c("LSC", "non LSC", "ESC")) + labs(x="tSNE 1", y="tSNE 2", title="")  + theme(plot.title = element_text(size=75, hjust = 0.5, face="bold"), legend.text = element_text(size=50))+ theme(aspect.ratio = 1) + theme(legend.title=element_blank())
dev.off()
# Shown in Extended Data Figure 3f 
