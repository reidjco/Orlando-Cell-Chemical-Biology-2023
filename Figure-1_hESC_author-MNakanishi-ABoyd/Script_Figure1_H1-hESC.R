################################################################################
# Orlando et al. Nature 2021
# R script used for scRNAseq H1 hESC data analysis shown in Extended Data Figure 1b-e,o 
# written by: Dr. Mio Nakanishi
# Date: October 3, 2019
################################################################################

# The count matrix file is publicly available on NCBI, GEO accession: GSE126022 
# Originally published: Nakanishi et al. Cell 2019. PMID: 30982595

################################################################################

# load libraries
library(DropletUtils)
library(scater)
library(scran)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(pheatmap)
library(gelnet)
library(dplyr)
library(gdata)
library(DT)
library(biomaRt)
library(synapseClient)

################################################################################

## NOTE: It is not necessary to re-derive the H1 count matrix, as it is publicly available on NCBI, GEO accession: GSE126022,
## and will be loaded in the next section. This section documents the separation of H1 (male) from H9 (female) hESCs

## LOADING THE COUNT MATRIX AND SEPARATION OF BULK H1 CELLS
## Load the count matrix
sce <- read10xCounts("~/folder_name", col.names = TRUE)
dim(sce)

BulkSex <- read.delim(file="sex_annotation_Bulk_H1_H9.txt", header=TRUE)
NCADposSex <- read.delim(file="sex_annotation_NCADpos_H1_H9.txt", header=TRUE)

BulkSex$Cell_Barcode <- paste(sep="",BulkSex$Cell_Barcode,"-1")
NCADposSex$Cell_Barcode <- paste(sep="",NCADposSex$Cell_Barcode,"-2")
MergedSex <- rbind(BulkSex, NCADposSex)

SexAnnotation<- MergedSex$Sex[match(colnames(sce), MergedSex$Cell_Barcode)]
sce$Sex <- SexAnnotation

##Separate H1(male) and H9(female). 
H1 <- filter(sce, sce$Sex == "male")
dim(H1)

## ESSENTIAL NOTE
## All analyses below will be performed on H1.
## H1 object will be renamed as 'sce'.
## Original H1 object will be removed to reduce the filesize.

sce <- H1
cellcodes <- as.data.frame(sce$Barcode)
colnames(cellcodes) <- "barcodes"
rownames(cellcodes) <- cellcodes$barcodes
cellcodes$libcodes <- as.factor(gsub(pattern=".+-", replacement="", cellcodes$barcode))
sce$grouping <- cellcodes$libcodes
table(sce$grouping)
sce <- filter(sce, sce$grouping == 1)
dim(sce)

################################################################################

## LOAD COUNT MATRIX FROM GEO AND CREATE STEM CELL EXPERIMENT

# set working directory
# download H1 count matrix from on NCBI, GEO accession: GSE126022

# load data
H1_GEO <- read.table("GSE126022_bulk_NCADpos_H1_counts.txt", header=TRUE,sep="\t")
# remove any potentially duplicated rows
H1_GEO <- H1_GEO[!duplicated(H1_GEO$X), ]
rownames(H1_GEO) <- H1_GEO[,1]
H1_GEO[,1] <- NULL # remove column X
H1_GEO[,1] <- NULL # remove column X.1

# set-up stem cell experiment (sce)
all.counts <- H1_GEO
all.counts <- as.matrix(all.counts)
cell_labels <- colnames(H1_GEO)
genes <- rownames(H1_GEO)

sce <- SingleCellExperiment(list(counts=all.counts),
                            colData=DataFrame(cell_labels),
                            rowData=DataFrame(genes),
                            metadata=list(study="H1_GEO"))

################################################################################

## QUALITY CONTROL OF THE CELLS
## Making mitochondrial gene list
is.mito <- grepl("^MT-", rowData(sce)$Symbol)
summary(is.mito)

## Generating QC Metrics
sce <- calculateQCMetrics(sce, feature_controls=list(Mt=is.mito)) 
par(mfrow=c(2,2), mar=c(5.1, 4.1, 0.1, 0.1))
hist(sce$total_counts/1e3, xlab="Library sizes (thousands)", main="", 
     breaks=20, col="grey80", ylab="Number of cells")
hist(sce$total_features_by_counts, xlab="Number of expressed genes", main="", 
     breaks=20, col="grey80", ylab="Number of cells")
hist(sce$pct_counts_Mt, xlab="Mitochondrial proportion (%)", 
     ylab="Number of cells", breaks=20, main="", col="grey80")

libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE)
feature.drop <- isOutlier(sce$total_features_by_counts, nmads=3, type="lower", log=TRUE)
mito.drop <- isOutlier(sce$pct_counts_Mt, nmads=3, type="higher", log=TRUE)
sce <- sce[,!(libsize.drop | feature.drop | mito.drop)]
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), ByMito=sum(mito.drop), Remaining=ncol(sce))

## Top 50 most abundant transcripts
fontsize <- theme(axis.text=element_text(size=12), axis.title=element_text(size=16))
plotQC(sce, type = "highest-expression", n=50) + fontsize

## Average counts for all genes
ave.counts <- calcAverage(sce, use_size_factors=FALSE)
hist(log10(ave.counts), breaks=100, main="", col="grey",
     xlab=expression(Log[10]~"average count"))

## Save the average counts for later use.
## Genes not expressed in any cells are deleted.
rowData(sce)$ave.count <- ave.counts
to.keep <- ave.counts > 0
sce <- sce[to.keep,]
summary(to.keep)

## Normalization of cell-specific bias
clusters <- quickCluster(sce, min.mean=0.1, method="igraph")
sce <- computeSumFactors(sce, cluster=clusters, min.mean=0.1)
summary(sizeFactors(sce))

plot(sizeFactors(sce), sce$total_counts/1e3, log="xy",
     ylab="Library size (thousands)", xlab="Size factor")
sce <- normalize(sce)

## Removing technical noise
var.fit.nospike <- trendVar(sce, parametric=TRUE, 
                            use.spikes=FALSE, loess.args=list(span=0.2))
var.out.nospike <- decomposeVar(sce, var.fit.nospike)

plot(var.out.nospike$mean, var.out.nospike$total, pch=16, cex=0.6, 
     xlab="Mean log-expression", ylab="Variance of log-expression")
curve(var.fit.nospike$trend(x), col="dodgerblue", lwd=2, add=TRUE)
points(var.out.nospike$mean[cur.spike], var.out.nospike$total[cur.spike], col="red", pch=16)

chosen.genes <- order(var.out.nospike$bio, decreasing=TRUE)[1:10]
plotExpression(sce, rownames(var.out.nospike)[chosen.genes], 
               alpha=0.05, jitter="jitter") + fontsize

sce <- denoisePCA(sce, technical=var.fit.nospike$trend, approximate=TRUE)
ncol(reducedDim(sce, "PCA"))

################################################################################

## CLUSTERING AND RE_SORTING

## Replace the rownames with gene symbols
rowData(sce)$ENSEMBL <- rownames(sce)
rownames(sce) <- rowData(sce)$Symbol

snn.gr <- buildSNNGraph(sce, use.dimred="PCA")
cluster.out <- igraph::cluster_walktrap(snn.gr)
my.clusters.igraph <- cluster.out$membership
table(my.clusters.igraph)

rm(cluster.out)
sce$cluster_igraph <- factor(my.clusters.igraph)
plotTSNE(sce, colour_by="cluster_igraph") + fontsize
 
##Draw PCA plots to sort the clusters.
plotReducedDim(sce, use_dimred="PCA", colour_by="cluster_igraph") + fontsize

## Draw modularity-heatmaps for the 2 methods.
mod.out.igraph <- clusterModularity(snn.gr, my.clusters.igraph, get.values=TRUE)
ratio.igraph <- log10(mod.out.igraph$observed/mod.out.igraph$expected + 1)
pheatmap(ratio.igraph, cluster_rows=FALSE, cluster_cols=FALSE, 
         color=colorRampPalette(c("white", "blue4"))(100))

## Reorder the clusters
cluster_igraph_reordered <- my.clusters.igraph
cluster_igraph_reordered <- paste(sep="",cluster_igraph_reordered,"a")
cluster_igraph_reordered <- gsub("10a", "1", cluster_igraph_reordered)
cluster_igraph_reordered <- gsub("7a", "2", cluster_igraph_reordered)
cluster_igraph_reordered <- gsub("13a", "8", cluster_igraph_reordered)
cluster_igraph_reordered <- gsub("3a", "3", cluster_igraph_reordered)
cluster_igraph_reordered <- gsub("12a", "4", cluster_igraph_reordered)
cluster_igraph_reordered <- gsub("6a", "5", cluster_igraph_reordered)
cluster_igraph_reordered <- gsub("2a", "6", cluster_igraph_reordered)
cluster_igraph_reordered <- gsub("5a", "7", cluster_igraph_reordered)

cluster_igraph_reordered <- gsub("4a", "9", cluster_igraph_reordered)
cluster_igraph_reordered <- gsub("11a", "10", cluster_igraph_reordered)
cluster_igraph_reordered <- gsub("1a", "11", cluster_igraph_reordered)
cluster_igraph_reordered <- gsub("8a", "12", cluster_igraph_reordered)
cluster_igraph_reordered <- gsub("9a", "13", cluster_igraph_reordered)
cluster_igraph_reordered <- as.numeric(cluster_igraph_reordered)
sce$cluster_igraph_reordered <- as.factor(cluster_igraph_reordered)

temp <- as.data.frame(reducedDim(sce, "TSNE"))
temp$col <- cluster_igraph_reordered
png("BulkH1_tSNE_reordered_igraph.png", width = 1000, height = 1000)
ggplot(temp, aes(x=V1, y=V2)) + geom_point(size=3, aes(color=factor(col))) + nogrids + theme(legend.position="none")
dev.off()
# Shown in Extended Data Figure 1j

################################################################################

## EXTENDED DATA FIGURE 1B-E: Gifford et al. Pluripotency and Lineage Scores

# Pluri score
pluri.markers <- scan("PluriMarkers_Meissner.txt", what="", sep="\t")
pluri.markers <- intersect(rownames(sce), pluri.markers)
pluri.vals <- as.matrix(logcounts(sce)[pluri.markers,,drop=FALSE])
all.vals <- as.matrix(logcounts(sce))

pluri.sig <- colMeans(pluri.vals)/nrow(pluri.vals) - colMeans(all.vals)/nrow(all.vals)
pluri.sig <- as.matrix(pluri.sig)
pluri.z <- pluri.sig/sd(pluri.sig) - mean(pluri.sig)/sd(pluri.sig)
colData(sce)$pluri.Z <- pluri.z

# Lineage scores
ecto.markers <- scan("EctoMarkers_Meissner.txt", what="", sep="\t")
ecto.markers <- intersect(rownames(sce), ecto.markers)
ecto.vals <- as.matrix(logcounts(sce)[ecto.markers,,drop=FALSE])
all.vals <- as.matrix(logcounts(sce))

ecto.sig <- colMeans(ecto.vals)/nrow(ecto.vals) - colMeans(all.vals)/nrow(all.vals)
ecto.sig <- as.matrix(ecto.sig)
ecto.z <- ecto.sig/sd(ecto.sig) - mean(ecto.sig)/sd(ecto.sig)
colData(sce)$ecto.Z <- ecto.z

meso.markers <- scan("MesoMarkers_Meissner.txt", what="", sep="\t")
meso.markers <- intersect(rownames(sce), meso.markers)
meso.vals <- as.matrix(logcounts(sce)[meso.markers,,drop=FALSE])

meso.sig <- colMeans(meso.vals)/nrow(meso.vals) - colMeans(all.vals)/nrow(all.vals)
meso.sig <- as.matrix(meso.sig)
meso.z <- meso.sig/sd(meso.sig) - mean(meso.sig)/sd(meso.sig)
colData(sce)$meso.Z <- meso.z

endo.markers <- scan("EndoMarkers_Meissner.txt", what="", sep="\t")
endo.markers <- intersect(rownames(sce), endo.markers)
endo.vals <- as.matrix(logcounts(sce)[endo.markers,,drop=FALSE])

endo.sig <- colMeans(endo.vals)/nrow(endo.vals) - colMeans(all.vals)/nrow(all.vals)
endo.sig <- as.matrix(endo.sig)
endo.z <- endo.sig/sd(endo.sig) - mean(endo.sig)/sd(endo.sig)
colData(sce)$endo.Z <- endo.z

## t-SNE PLOTS
collist<-c("cluster_igraph_reordered", "ecto.Z", "meso.Z", "endo.Z", "pluri.Z")
temp <- colData(sce)[,collist]

temp$ecto.Z <- as.numeric(temp$ecto.Z)
temp$endo.Z <- as.numeric(temp$endo.Z)
temp$meso.Z <- as.numeric(temp$meso.Z)
temp$pluri.Z <- as.numeric(temp$pluri.Z)
temp <- as.data.frame(temp)

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

temp <- as.data.frame(reducedDim(sce, "TSNE"))
temp$col <- as.numeric(sce$pluri.Z)
png("BulkH9_histogram_pluriZ3.png", width = 1000, height = 1000)
ggplot(temp, aes(x=col)) + geom_density(color="firebrick3", fill="firebrick1", lwd = 2.5) + theme_classic(base_size = 36) + labs(x="PLURI score", y="Cell number", title="PLURI score") + theme(plot.title = element_text(size=50, hjust = 0.5, face="bold")) + theme(aspect.ratio = 1) +geom_vline(xintercept = 0.8, size = 2, linetype = "longdash", alpha = 0.75)
dev.off()

nogrids <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))

# t-SNE coloured by Gifford Pluri gene signature, Extended Data Figure 1b
temp <- as.data.frame(reducedDim(sce, "TSNE"))
temp$col <- pluri.z
png("BulkH9_tSNE_pluriZ2_wlegend.png", width = 1000, height = 1000)
ggplot(temp, aes(x=V1, y=V2)) + geom_point(size=4, aes(color=col)) + theme_classic(base_size = 50) + labs(x="tSNE 1", y="tSNE 2", title="PLURI score")  + scale_colour_gradient2(low = "blue3", mid = "snow2", high = "red2", midpoint = mean(temp$col), space = "Lab", guide = "colourbar", aesthetics = "colour") + theme(plot.title = element_text(size=75, hjust = 0.5, face="bold"))+ theme(aspect.ratio = 1) + theme(legend.title=element_blank()) + theme(legend.text = element_text(size = 20))
dev.off()

# t-SNE coloured by Gifford Meso gene signature, Extended Data Figure 1c
temp <- as.data.frame(reducedDim(sce, "TSNE"))
temp$col <- meso.z
png("BulkH9_tSNE_mesoZ2.png", width = 1000, height = 1000)
ggplot(temp, aes(x=V1, y=V2)) + geom_point(size=4, aes(color=col)) + theme_classic(base_size = 50) + labs(x="tSNE 1", y="tSNE 2", title="MESO score")  + scale_colour_gradient2(low = "blue3", mid = "snow2", high = "red2", midpoint = mean(temp$col), space = "Lab", guide = "colourbar", aesthetics = "colour") + theme(plot.title = element_text(size=75, hjust = 0.5, face="bold"))+ theme(aspect.ratio = 1) + theme(legend.position="none") 
dev.off()

# t-SNE coloured by Gifford Ecto gene signature, Extended Data Figure 1d
temp <- as.data.frame(reducedDim(sce, "TSNE"))
temp$col <- ecto.z
png("BulkH9_tSNE_ectoZ2.png", width = 1000, height = 1000)
ggplot(temp, aes(x=V1, y=V2)) + geom_point(size=4, aes(color=col)) + theme_classic(base_size = 50) + labs(x="tSNE 1", y="tSNE 2", title="ECTO score")  + scale_colour_gradient2(low = "blue3", mid = "snow2", high = "red2", midpoint = mean(temp$col), space = "Lab", guide = "colourbar", aesthetics = "colour") + theme(plot.title = element_text(size=75, hjust = 0.5, face="bold"))+ theme(aspect.ratio = 1) + theme(legend.position="none") 
dev.off()

# t-SNE coloured by Gifford Endo gene signature, Extended Data Figure 1e
temp <- as.data.frame(reducedDim(sce, "TSNE"))
temp$col <- endo.z
png("BulkH9_tSNE_endoZ2.png", width = 1000, height = 1000)
ggplot(temp, aes(x=V1, y=V2)) + geom_point(size=4, aes(color=col)) + theme_classic(base_size = 50) + labs(x="tSNE 1", y="tSNE 2", title="ENDO score")  + scale_colour_gradient2(low = "blue3", mid = "snow2", high = "red2", midpoint = mean(temp$col), space = "Lab", guide = "colourbar", aesthetics = "colour") + theme(plot.title = element_text(size=75, hjust = 0.5, face="bold"))+ theme(aspect.ratio = 1) + theme(legend.position="none") 
dev.off()

## VIOLIN PLOTS
collist<-c("cluster_igraph_reordered", "ecto.Z", "meso.Z", "endo.Z", "pluri.Z")
temp <- colData(sce)[,collist]

temp$ecto.Z <- as.numeric(temp$ecto.Z)
temp$endo.Z <- as.numeric(temp$endo.Z)
temp$meso.Z <- as.numeric(temp$meso.Z)
temp$pluri.Z <- as.numeric(temp$pluri.Z)
temp <- as.data.frame(temp)

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# violin plot scored by Gifford Pluri gene signature, Extended Data Figure 1b
png("violin_pluriZ_reordered_Igraph.png", width = 1000, height = 1000)
ggplot(temp, aes(x=cluster_igraph_reordered, y=pluri.Z, fill=cluster_igraph_reordered)) + nogrids + geom_violin() + stat_summary(fun.data=data_summary) + theme(legend.position="none")
dev.off()

# violin plot scored by Gifford Meso gene signature, Extended Data Figure 1c
png("violin_mesoZ_reordered_Igraph.png", width = 1000, height = 1000)
ggplot(temp, aes(x=cluster_igraph_reordered, y=meso.Z, fill=cluster_igraph_reordered)) + nogrids + geom_violin() + stat_summary(fun.data=data_summary) + theme(legend.position="none")
dev.off()

# violin plot scored by Gifford Ecto gene signature, Extended Data Figure 1d
png("violin_ectoZ_reordered_Igraph.png", width = 1000, height = 1000)
ggplot(temp, aes(x=cluster_igraph_reordered, y=ecto.Z, fill=cluster_igraph_reordered)) + nogrids + geom_violin() + stat_summary(fun.data=data_summary) + theme(legend.position="none")
dev.off()

# violin plot scored by Gifford Endo gene signature, Extended Data Figure 1e
png("violin_endoZ_reordered_Igraph.png", width = 1000, height = 1000)
ggplot(temp, aes(x=cluster_igraph_reordered, y=endo.Z, fill=cluster_igraph_reordered)) + nogrids + geom_violin() + stat_summary(fun.data=data_summary) + theme(legend.position="none")
dev.off()

################################################################################

## EXTENDED DATA FIGURE 1o: Stem Cell Index Score

temp <- logcounts(sce)
temp <- as.matrix(temp)
temp <- temp[!duplicated(row.names(temp)), ]
write.table(temp, file = "BulkH1_logcounts.txt", sep="\t", row.names=TRUE, col.names=TRUE)
rm(temp)

w <- read.delim("pcbc-stemsig.tsv", header=FALSE, row.names=1) %>% as.matrix() %>% drop()

## Reduces HUGO|POSITION gene IDs to just HUGO
f <- function( v ) unlist( lapply( strsplit( v, "\\|" ), "[[", 1 ) )

X <- read.delim("BulkH1_logcounts.txt", as.is=TRUE, check.names=FALSE, sep="\t") %>%  ## Read the raw values
  filter( !grepl( "\\?", gene_id ) ) %>%      ## Drop genes with no mapping to HUGO
  mutate( gene_id = f( gene_id ) ) %>%        ## Clip gene ids to HUGO
  filter( gene_id %in% names(w) )         ## Reduce to the signature's gene set

## SLC35E2 has multiple entries with the same HUGO id
## Keep the first entry only
j <- grep( "SLC35E2", X[,1] )
if( length(j) > 1 )
  X <- X[-j[-1],]

## Convert to a matrix
rownames(X) <- NULL
X <- X %>% tibble::column_to_rownames( "gene_id" ) %>% as.matrix()

## Reduce the signature to the common set of genes
stopifnot( all( rownames(X) %in% names(w) ) )
w <- w[ rownames(X) ]

####### Score via Spearman correlation
s <- apply( X, 2, function(z) {cor( z, w, method = "sp", use = "complete.obs" )} )

## Scale the scores to be between 0 and 1
s <- s - min(s)
s <- s / max(s)

write.table(cbind(s), file = "BulkH1_scRNA-seq_StemScore.tsv", sep = "\t", quote = FALSE, col.names = FALSE)

mRNAsi <- read.delim("BulkH1_scRNA-seq_StemScore.tsv", header=FALSE, row.names=1)
cellID <- colnames(sce)
mRNAsi <- mRNAsi$V2[match(cellID, rownames(mRNAsi))]
mRNAsi <- t(mRNAsi)
colData(sce)$mRNAsi <- mRNAsi



png("BulkH1_violin_mRNAsi_reordered_Igraph.png", width = 1000, height = 1000)
ggplot(temp, aes(x=cluster_igraph_reordered, y=mRNAsi, fill=cluster_igraph_reordered)) + nogrids + geom_violin() + stat_summary(fun.data=data_summary) + theme(legend.position="none")
dev.off()

################################################################################

## statistics for Extended Data Figure 1b-e
difscores <- colData(sce)[,c("cluster_igraph_reordered", "ecto.Z", "meso.Z", "endo.Z", "pluri.Z")]
difscores <- as.data.frame(difscores)

pluri.model <- aov(pluri.Z~cluster_igraph_reordered, data=difscores)
summary(pluri.model) 
TukeyHSD(pluri.model, conf.level = 0.99)