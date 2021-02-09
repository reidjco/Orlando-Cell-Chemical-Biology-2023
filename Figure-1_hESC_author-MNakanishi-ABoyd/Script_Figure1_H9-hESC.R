################################################################################
# Orlando et al. Nature 2021
# R script used for scRNAseq H9 hESC data analysis 
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

## LOADING THE COUNT MATRIX AND SEPARATION OF BULK H9 CELLS.
## Load the count matrix
sce <- read10xCounts("~/folder_name", col.names = TRUE)
dim(sce)

## Load sex-annotation-file. "Cell Barcode" header was changed to "Cell_Barcode" using Text Editor.
BulkSex <- read.delim(file="sex_annotation_Bulk_H1_H9.txt", header=TRUE)
NCADposSex <- read.delim(file="sex_annotation_NCADpos_H1_H9.txt", header=TRUE)

## In the merged count matrix (and sce object), "-1" and "-2" were attached to the end of cell barcode
## (colnames of sce) for the Bulk and NCADpos libraries, respectively. Modify the cell barcodes in the 
## retrieved sex-annotation objects accordingly. 

BulkSex$Cell_Barcode <- paste(sep="",BulkSex$Cell_Barcode,"-1")
NCADposSex$Cell_Barcode <- paste(sep="",NCADposSex$Cell_Barcode,"-2")
MergedSex <- rbind(BulkSex, NCADposSex)

## Sort the merged sex-annotation by the 
## order of colnames(sce) (cell barcodes in the matrix) and put them sce$Sex. 
SexAnnotation<- MergedSex$Sex[match(colnames(sce), MergedSex$Cell_Barcode)]
sce$Sex <- SexAnnotation

## Separate H1(male) and H9(female). 
H9 <- filter(sce, sce$Sex == "female")
dim(H9)

## ESSENTIAL NOTE
## All analyses below will be performed on H9.
## H9 object will be renamed as 'sce'.
## Original H9 object will be removed to reduce the filesize.
sce <- H9

## Furthermore, separate 2 libraries (grouping (1: Bulk, 2: NCADpos)) and remove NCADpos.
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
# download H9 count matrix from on NCBI, GEO accession: GSE126022

# load data
H9_GEO <- read.table("GSE126022_bulk_NCADpos_H9_counts.txt", header=TRUE,sep="\t")
# remove any potentially duplicated rows
H9_GEO <- H9_GEO[!duplicated(H9_GEO$X), ]
rownames(H9_GEO) <- H9_GEO[,1]
H9_GEO[,1] <- NULL # remove column X
H9_GEO[,1] <- NULL # remove column X.1

# set-up stem cell experiment (sce)
all.counts <- H9_GEO
all.counts <- as.matrix(all.counts)
cell_labels <- colnames(H9_GEO)
genes <- rownames(H9_GEO)

sce <- SingleCellExperiment(list(counts=all.counts),
                            colData=DataFrame(cell_labels),
                            rowData=DataFrame(genes),
                            metadata=list(study="H9_GEO"))

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
hist(sce$total_features, xlab="Number of expressed genes", main="", 
     breaks=20, col="grey80", ylab="Number of cells")
hist(sce$pct_counts_Mt, xlab="Mitochondrial proportion (%)", 
     ylab="Number of cells", breaks=20, main="", col="grey80")

libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE)
feature.drop <- isOutlier(sce$total_features, nmads=3, type="lower", log=TRUE)
mito.drop <- isOutlier(sce$pct_counts_Mt, nmads=3, type="higher", log=TRUE)
sce <- sce[,!(libsize.drop | feature.drop | mito.drop)]
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), ByMito=sum(mito.drop), Remaining=ncol(sce))

## Cell cycle analysis
hs.pairs <- readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
assignments <- cyclone(sce, hs.pairs, gene.names=rownames(sce))
sce$phases <- assignments$phases
table(sce$phases)
plot(assignments$score$G1, assignments$score$G2M, xlab="G1 score", ylab="G2/M score", pch=16)

## EXAMININE GENE-LEVEL METRICS
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

## NORMALIZATION OF CELL-SPECIFIC BIASES
## computeSumFactors
clusters <- quickCluster(sce, min.mean=0.1, method="igraph")
sce <- computeSumFactors(sce, cluster=clusters, min.mean=0.1)
summary(sizeFactors(sce))

plot(sizeFactors(sce), sce$total_counts/1e3, log="xy",
     ylab="Library size (thousands)", xlab="Size factor")

sce <- normalize(sce)

## MODELLING AND REMOVING TECHNICAL NOISE
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

## CLUSTERING 
## Replace the rownames with gene symbols
rowData(sce)$ENSEMBL <- rownames(sce)
rownames(sce) <- rowData(sce)$Symbol

## Test clustering method iGraph.
## These cluster#s were re-sorted in following section
snn.gr <- buildSNNGraph(sce, use.dimred="PCA")
cluster.out <- igraph::cluster_walktrap(snn.gr)
my.clusters.igraph <- cluster.out$membership
table(my.clusters.igraph)

sce$cluster_igraph <- factor(my.clusters.igraph)
plotTSNE(sce, colour_by="cluster_igraph") + fontsize

## Draw PCA plots to sort the clusters.
plotReducedDim(sce, use_dimred="PCA", colour_by="cluster_igraph") + fontsize

## Draw modularity-heatmaps for the 2 methods.
mod.out.igraph <- clusterModularity(snn.gr, my.clusters.igraph, get.values=TRUE)
ratio.igraph <- log10(mod.out.igraph$observed/mod.out.igraph$expected + 1)
pheatmap(ratio.igraph, cluster_rows=FALSE, cluster_cols=FALSE, 
         color=colorRampPalette(c("white", "blue4"))(100))

################################################################################

## RE-SORTING THE CLUSTERS
## Given the gradation from primitive to differentiated clusters
## appeared in the t-SNE colored by the Stem Cell Index, the cluster_igraph were re-sorted.
cluster_igraph_reordered <- my.clusters.igraph
cluster_igraph_reordered <- paste(sep="",cluster_igraph_reordered,"a")
cluster_igraph_reordered <- gsub("10a", "1", cluster_igraph_reordered)
cluster_igraph_reordered <- gsub("7a", "2", cluster_igraph_reordered)
cluster_igraph_reordered <- gsub("6a", "3", cluster_igraph_reordered)
cluster_igraph_reordered <- gsub("4a", "4", cluster_igraph_reordered)
cluster_igraph_reordered <- gsub("5a", "5", cluster_igraph_reordered)
cluster_igraph_reordered <- gsub("3a", "6", cluster_igraph_reordered)
cluster_igraph_reordered <- gsub("8a", "7", cluster_igraph_reordered)
cluster_igraph_reordered <- gsub("1a", "8", cluster_igraph_reordered)
cluster_igraph_reordered <- gsub("9a", "9", cluster_igraph_reordered)
cluster_igraph_reordered <- gsub("2a", "10", cluster_igraph_reordered)
cluster_igraph_reordered <- as.numeric(cluster_igraph_reordered)
sce$cluster_igraph_reordered <- as.factor(cluster_igraph_reordered)

################################################################################

## FIGURE 1A,G: t-SNE plots of H9 hESC reordered clusters
temp$col <- cluster_igraph_reordered
png("tSNE_reordered_igraph.png", width = 1000, height = 1000)
ggplot(temp, aes(x=V1, y=V2)) + geom_point(size=3, aes(color=factor(col))) + nogrids + theme(legend.position="none")
dev.off()

################################################################################

## FIGURE 1B,D-F: Gifford et al. Pluripotency and Lineage Scores

# Pluri score
pluri.markers <- scan("List_Gifford_PluriMarkers.txt", what="", sep="\t")
pluri.markers <- intersect(rownames(sce), pluri.markers)
pluri.vals <- as.matrix(logcounts(sce)[pluri.markers,,drop=FALSE])
all.vals <- as.matrix(logcounts(sce))

pluri.sig <- colMeans(pluri.vals)/nrow(pluri.vals) - colMeans(all.vals)/nrow(all.vals)
pluri.sig <- as.matrix(pluri.sig)
pluri.z <- pluri.sig/sd(pluri.sig) - mean(pluri.sig)/sd(pluri.sig)
colData(sce)$pluri.Z <- pluri.z

# Lineage scores
ecto.markers <- scan("List_Gifford_EctoMarkers.txt", what="", sep="\t")
ecto.markers <- intersect(rownames(sce), ecto.markers)
ecto.vals <- as.matrix(logcounts(sce)[ecto.markers,,drop=FALSE])
all.vals <- as.matrix(logcounts(sce))

ecto.sig <- colMeans(ecto.vals)/nrow(ecto.vals) - colMeans(all.vals)/nrow(all.vals)
ecto.sig <- as.matrix(ecto.sig)
ecto.z <- ecto.sig/sd(ecto.sig) - mean(ecto.sig)/sd(ecto.sig)
colData(sce)$ecto.Z <- ecto.z

meso.markers <- scan("List_Gifford_MesoMarkers.txt", what="", sep="\t")
meso.markers <- intersect(rownames(sce), meso.markers)
meso.vals <- as.matrix(logcounts(sce)[meso.markers,,drop=FALSE])

meso.sig <- colMeans(meso.vals)/nrow(meso.vals) - colMeans(all.vals)/nrow(all.vals)
meso.sig <- as.matrix(meso.sig)
meso.z <- meso.sig/sd(meso.sig) - mean(meso.sig)/sd(meso.sig)
colData(sce)$meso.Z <- meso.z

endo.markers <- scan("List_Gifford_EndoMarkers.txt", what="", sep="\t")
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

# t-SNE coloured by Gifford Pluri gene signature, Figure 1b
temp <- as.data.frame(reducedDim(sce, "TSNE"))
temp$col <- pluri.z
png("BulkH9_tSNE_pluriZ2_wlegend.png", width = 1000, height = 1000)
ggplot(temp, aes(x=V1, y=V2)) + geom_point(size=4, aes(color=col)) + theme_classic(base_size = 50) + labs(x="tSNE 1", y="tSNE 2", title="PLURI score")  + scale_colour_gradient2(low = "blue3", mid = "snow2", high = "red2", midpoint = mean(temp$col), space = "Lab", guide = "colourbar", aesthetics = "colour") + theme(plot.title = element_text(size=75, hjust = 0.5, face="bold"))+ theme(aspect.ratio = 1) + theme(legend.title=element_blank()) + theme(legend.text = element_text(size = 20))
dev.off()

# t-SNE coloured by Gifford Meso gene signature, Figure 1d
temp <- as.data.frame(reducedDim(sce, "TSNE"))
temp$col <- meso.z
png("BulkH9_tSNE_mesoZ2.png", width = 1000, height = 1000)
ggplot(temp, aes(x=V1, y=V2)) + geom_point(size=4, aes(color=col)) + theme_classic(base_size = 50) + labs(x="tSNE 1", y="tSNE 2", title="MESO score")  + scale_colour_gradient2(low = "blue3", mid = "snow2", high = "red2", midpoint = mean(temp$col), space = "Lab", guide = "colourbar", aesthetics = "colour") + theme(plot.title = element_text(size=75, hjust = 0.5, face="bold"))+ theme(aspect.ratio = 1) + theme(legend.position="none") 
dev.off()

# t-SNE coloured by Gifford Ecto gene signature, Figure 1e
temp <- as.data.frame(reducedDim(sce, "TSNE"))
temp$col <- ecto.z
png("BulkH9_tSNE_ectoZ2.png", width = 1000, height = 1000)
ggplot(temp, aes(x=V1, y=V2)) + geom_point(size=4, aes(color=col)) + theme_classic(base_size = 50) + labs(x="tSNE 1", y="tSNE 2", title="ECTO score")  + scale_colour_gradient2(low = "blue3", mid = "snow2", high = "red2", midpoint = mean(temp$col), space = "Lab", guide = "colourbar", aesthetics = "colour") + theme(plot.title = element_text(size=75, hjust = 0.5, face="bold"))+ theme(aspect.ratio = 1) + theme(legend.position="none") 
dev.off()

# t-SNE coloured by Gifford Endo gene signature, Figure 1f
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

# violin plot scored by Gifford Pluri gene signature, Figure 1b
png("violin_pluriZ_reordered_Igraph.png", width = 1000, height = 1000)
ggplot(temp, aes(x=cluster_igraph_reordered, y=pluri.Z, fill=cluster_igraph_reordered)) + nogrids + geom_violin() + stat_summary(fun.data=data_summary) + theme(legend.position="none")
dev.off()

# violin plot scored by Gifford Meso gene signature, Figure 1d
png("violin_mesoZ_reordered_Igraph.png", width = 1000, height = 1000)
ggplot(temp, aes(x=cluster_igraph_reordered, y=meso.Z, fill=cluster_igraph_reordered)) + nogrids + geom_violin() + stat_summary(fun.data=data_summary) + theme(legend.position="none")
dev.off()

# violin plot scored by Gifford Ecto gene signature, Figure 1e
png("violin_ectoZ_reordered_Igraph.png", width = 1000, height = 1000)
ggplot(temp, aes(x=cluster_igraph_reordered, y=ecto.Z, fill=cluster_igraph_reordered)) + nogrids + geom_violin() + stat_summary(fun.data=data_summary) + theme(legend.position="none")
dev.off()

# violin plot scored by Gifford Endo gene signature, Figure 1f
png("violin_endoZ_reordered_Igraph.png", width = 1000, height = 1000)
ggplot(temp, aes(x=cluster_igraph_reordered, y=endo.Z, fill=cluster_igraph_reordered)) + nogrids + geom_violin() + stat_summary(fun.data=data_summary) + theme(legend.position="none")
dev.off()

################################################################################

## FIGURE 1C: OCT4 Expression

## t-SNE coloured by OCT4 (gene name is POU5F1) single gene's expression level, Figure 1c
temp <- as.data.frame(reducedDim(sce, "TSNE"))
temp$col <- logcounts(sce)["POU5F1",]
png("BulkH9_tSNE_OCT4.png", width = 1000, height = 1000)
ggplot(temp, aes(x=V1, y=V2)) + geom_point(size=3, aes(color=col)) + theme_classic(base_size = 40) + labs(x="tSNE 1", y="tSNE 2", title="OCT4")  + scale_colour_gradient(low = "grey85", high = "red", guide = "colourbar", aesthetics = "colour") + theme(plot.title = element_text(size=75, hjust = 0.5, face="bold"))+ theme(aspect.ratio = 1) + theme(legend.position="none") 
dev.off()

################################################################################

## EXTENDED DATA FIGURE 1A: Pluripotency Scores

# Sperger et al. PLURI score, Extended Data Figure 1a (left)
temp.markers <- scan("~/List_Sperger_PluriMarkers.txt", what="", sep="\t")
temp.markers <- intersect(rownames(sce), temp.markers)
temp.vals <- as.matrix(logcounts(sce)[temp.markers,,drop=FALSE])
temp.sig <- colMeans(temp.vals)/nrow(temp.vals) - colMeans(all.vals)/nrow(all.vals)
temp.sig <- as.matrix(temp.sig)
temp.z <- temp.sig/sd(temp.sig) - mean(temp.sig)/sd(temp.sig)
temp <- as.data.frame(temp.z)
temp$cluster_igraph_reordered <- as.factor(cluster_igraph_reordered)

png("violin_Sperger_cluster_igraph_reordered.png", width = 900, height = 1200)
ggplot(temp, aes(x=cluster_igraph_reordered, y=V1, fill=cluster_igraph_reordered)) + geom_violin(scale = "width", trim = FALSE, lwd = 1.5) + theme_classic(base_size = 36) + labs(x="Clusters", y="PLURI scores", title="PLURI score (Sperger et al.)") + theme(plot.title = element_text(size=58, hjust = 0.5, face="bold"), axis.text.x = element_text(face="bold", size=32)) + theme(aspect.ratio = 1.33) + theme(legend.position="none") 
dev.off()

# Hutchins et al. PLURI score using 'embryonic gene list', Extended Data Figure 1a (right)
all.vals <- as.matrix(logcounts(sce))

Embryonic.markers <- scan("List_Hutchins_PluriMarkers.txt", what="", sep="\t")
Embryonic.markers <- intersect(rownames(sce), Embryonic.markers)
Embryonic.vals <- as.matrix(logcounts(sce)[Embryonic.markers,,drop=FALSE])
Embryonic.sig <- colMeans(Embryonic.vals)/nrow(Embryonic.vals) - colMeans(all.vals)/nrow(all.vals)
Embryonic.sig <- as.matrix(Embryonic.sig)
Embryonic.z <- Embryonic.sig/sd(Embryonic.sig) - mean(Embryonic.sig)/sd(Embryonic.sig)
colData(sce)$Embryonic.Z <- Embryonic.z

png("violin_Hutchins_cluster_igraph_reordered.png", width = 900, height = 1200)
ggplot(temp, aes(x=cluster_igraph_reordered, y=Embryonic.z, fill=cluster_igraph_reordered)) + geom_violin(scale = "width", trim = FALSE, lwd = 1.5) + theme_classic(base_size = 36) + labs(x="Clusters", y="PLURI scores", title="PLURI score (Hutchins et al.)") + theme(plot.title = element_text(size=58, hjust = 0.5, face="bold"), axis.text.x = element_text(face="bold", size=32)) + theme(aspect.ratio = 1.33) + theme(legend.position="none") 
dev.off()

################################################################################

## EXTENDED DATA FIGURE 1F-H: Hutchins et al. Lineage Scores

## Retrieved gene symbols for the domain-specific genes in the Table S3: Markers_Hutchins.xlsx
## Save the gene symbols as tab-delimited .txt file, and store in working directory

Markers_Hutchins <- read.delim(file="Markers_Hutchins.txt", header=TRUE)
head(Markers_Hutchins)
##   Blood_mesoderm    Embryonic Endoderm    Germ_cells Mesoderm Neural_crest Neurectoderm Surface_ectoderm
## 1          Pydc3        Peg10    Peg10         Armc3    Peg10        Peg10       Begain            Peg10
## 2       BC051077        Cldn6    Cldn6       Col22a1   Adhfe1       C79798        Peg10            Bpifc
## 3         Rnf128 RP24-114A9.1   Col4a6 1700010I14Rik   Igfbp6    Hist1h2ac         Rgma             Rgma
## 4        Slc46a3      Trim43b   Rnf128        Tcp10b     Rgma         Rgma         Grm7           Igfbp7
## 5             Xk     BC051077   Igfbp7 4931429I11Rik   Rnf128      Gm26718         Fjx1             Sdpr
## 6          Cep72       Zfp534     Sdpr        Dnah7c    Fndc1         Fjx1         Gnaz              Dct

# Extract mouse genes
Embryonic.Hutchins <- Markers_Hutchins$Embryonic
Mesoderm.Hutchins <- Markers_Hutchins$Mesoderm
Endoderm.Hutchins <- Markers_Hutchins$Endoderm
Neurectoderm.Hutchins <- Markers_Hutchins$Neurectoderm

# Function to convert mouse to human genes:
convertMmtoHs <- function(x){
  require("biomaRt")
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  human.x <<- unique(genesV2[, 2])
}

# Call function to create: List_Hutchins_PluriMarkers.txt
convertMmtoHs(Embryonic.Hutchins)
human.Embryonic.Hutchins <- human.x
rm(human.x)
write.table(human.Embryonic.Hutchins, file = "List_Hutchins_PluriMarkers.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

name <- deparse(substitute(x))
name <- paste(sep="", "human_", name)
print(name) <<- human.x
dev.off()

# Call function to create: List_Hutchins_MesoMarkers.txt
convertMmtoHs(Mesoderm.Hutchins)
human.Mesoderm.Hutchins <- human.x
rm(human.x)
write.table(human.Mesoderm.Hutchins, file = "List_Hutchins_MesoMarkers.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

name <- deparse(substitute(x))
name <- paste(sep="", "human_", name)
print(name) <<- human.x
dev.off()

# Call function to create: List_Hutchins_EndoMarkers.txt
convertMmtoHs(Endoderm.Hutchins)
human.Endoderm.Hutchins <- human.x
rm(human.x)
write.table(human.Endoderm.Hutchins, file = "List_Hutchins_EndoMarkers.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

name <- deparse(substitute(x))
name <- paste(sep="", "human_", name)
print(name) <<- human.x
dev.off()

# Call function to create: List_Hutchins_EctoMarkers.txt
convertMmtoHs(Ectoderm.Hutchins)
human.Ectoderm.Hutchins <- human.x
rm(human.x)
write.table(human.Ectoderm.Hutchins, file = "List_Hutchins_EctoMarkers.txt", sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE)

name <- deparse(substitute(x))
name <- paste(sep="", "human_", name)
print(name) <<- human.x
dev.off()

## Score H9 scRNAseq using the Hutchins-derived human gene lists:
# Pluri score
pluri.markersHutchins <- scan("List_Hutchins_PluriMarkers.txt", what="", sep="\t")
pluri.markersHutchins <- intersect(rownames(sce), pluri.markersHutchins)
pluri.valsHutchins <- as.matrix(logcounts(sce)[pluri.markersHutchins,,drop=FALSE])
all.valsHutchins <- as.matrix(logcounts(sce))

pluri.sigHutchins <- colMeans(pluri.valsHutchins)/nrow(pluri.valsHutchins) - colMeans(all.valsHutchins)/nrow(all.valsHutchins)
pluri.sigHutchins <- as.matrix(pluri.sigHutchins)
pluri.zHutchins <- pluri.sigHutchins/sd(pluri.sigHutchins) - mean(pluri.sigHutchins)/sd(pluri.sigHutchins)
colData(sce)$pluri.zHutchins <- pluri.zHutchins

# Ecto score
ecto.markersHutchins <- scan("List_Hutchins_EctoMarkers.txt", what="", sep="\t")
ecto.markersHutchins <- intersect(rownames(sce), ecto.markersHutchins)
ecto.valsHutchins <- as.matrix(logcounts(sce)[ecto.markersHutchins,,drop=FALSE])
all.valsHutchins <- as.matrix(logcounts(sce))

ecto.sigHutchins <- colMeans(ecto.valsHutchins)/nrow(ecto.valsHutchins) - colMeans(all.valsHutchins)/nrow(all.valsHutchins)
ecto.sigHutchins <- as.matrix(ecto.sigHutchins)
ecto.zHutchins <- ecto.sigHutchins/sd(ecto.sigHutchins) - mean(ecto.sigHutchins)/sd(ecto.sigHutchins)
colData(sce)$ecto.zHutchins <- ecto.zHutchins

# Meso score
meso.markersHutchins <- scan("List_Hutchins_MesoMarkers.txt", what="", sep="\t")
meso.markersHutchins <- intersect(rownames(sce), meso.markersHutchins)
meso.valsHutchins <- as.matrix(logcounts(sce)[meso.markersHutchins,,drop=FALSE])

meso.sigHutchins <- colMeans(meso.valsHutchins)/nrow(meso.valsHutchins) - colMeans(all.valsHutchins)/nrow(all.valsHutchins)
meso.sigHutchins <- as.matrix(meso.sigHutchins)
meso.zHutchins <- meso.sigHutchins/sd(meso.sigHutchins) - mean(meso.sigHutchins)/sd(meso.sigHutchins)
colData(sce)$meso.zHutchins <- meso.zHutchins

# Endo score
endo.markersHutchins <- scan("List_Hutchins_EndoMarkers.txt", what="", sep="\t")
endo.markersHutchins <- intersect(rownames(sce), endo.markersHutchins)
endo.valsHutchins <- as.matrix(logcounts(sce)[endo.markersHutchins,,drop=FALSE])

endo.sigHutchins <- colMeans(endo.valsHutchins)/nrow(endo.valsHutchins) - colMeans(all.valsHutchins)/nrow(all.valsHutchins)
endo.sigHutchins <- as.matrix(endo.sigHutchins)
endo.zHutchins <- endo.sigHutchins/sd(endo.sigHutchins) - mean(endo.sigHutchins)/sd(endo.sigHutchins)
colData(sce)$endo.zHutchins <- endo.zHutchins

# t-SNE coloured by Hutchins Meso gene signature, Extended Data Figure 1f
temp <- as.data.frame(reducedDim(sce, "TSNE"))
temp$col <- meso.zHutchins
png("tSNE_Hutchins_MesoScore.png", width = 1000, height = 1000)
ggplot(temp, aes(x=V1, y=V2)) + geom_point(size=3, aes(color=col)) + nogrids + scale_colour_gradient2(low = "blue", mid = "yellow", high = "red", midpoint = max(temp$col)/2, space = "Lab", guide = "colourbar", aesthetics = "colour")  + theme(legend.position="none")
dev.off()

# t-SNE coloured by Hutchins Ecto gene signature, Extended Data Figure 1g
temp <- as.data.frame(reducedDim(sce, "TSNE"))
temp$col <- ecto.zHutchins
png("tSNE_Hutchins_EctoScore.png", width = 1000, height = 1000)
ggplot(temp, aes(x=V1, y=V2)) + geom_point(size=3, aes(color=col)) + nogrids + scale_colour_gradient2(low = "blue", mid = "yellow", high = "red", midpoint = max(temp$col)/2, space = "Lab", guide = "colourbar", aesthetics = "colour") + theme(legend.position="none")
dev.off()

# t-SNE coloured by Hutchins Endo gene signature, Extended Data Figure 1h
temp <- as.data.frame(reducedDim(sce, "TSNE"))
temp$col <- endo.zHutchins
png("tSNE_Hutchins_EndoScore.png", width = 1000, height = 1000)
ggplot(temp, aes(x=V1, y=V2)) + geom_point(size=3, aes(color=col)) + nogrids + scale_colour_gradient2(low = "blue", mid = "yellow", high = "red", midpoint = max(temp$col)/2, space = "Lab", guide = "colourbar", aesthetics = "colour")  + theme(legend.position="none")
dev.off()

################################################################################

## EXTENDED DATA FIGURE 1N: CALCULATE STEM CELL INDEX 

temp <- logcounts(sce)
temp <- as.matrix(temp)
temp <- temp[!duplicated(row.names(temp)), ]
write.table(temp, file = "BulkH9_logcounts.txt", sep="\t", row.names=TRUE, col.names=TRUE)
## Add "gene_id(tab)" at the head of the "BulkH9_logcounts.txt"  

w <- read.delim("pcbc-stemsig.tsv", header=FALSE, row.names=1) %>% as.matrix() %>% drop()

## Reduces HUGO|POSITION gene IDs to just HUGO
f <- function( v ) unlist( lapply( strsplit( v, "\\|" ), "[[", 1 ) )

X <- read.delim("BulkH9_logcounts.txt", as.is=TRUE, check.names=FALSE, sep="\t") %>%  ## Read the raw values
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

write.table(cbind(s), file = "BulkH9_scRNA-seq_StemScore.tsv", sep = "\t", quote = FALSE, col.names = FALSE)

# t-SNE plot, shown in Extended Data Figure 1n (left)
mRNAsi <- read.delim("BulkH9_scRNA-seq_StemScore.tsv", header=FALSE, row.names=1)
cellID <- colnames(sce)
mRNAsi <- mRNAsi$V2[match(cellID, rownames(mRNAsi))]
mRNAsi <- t(mRNAsi)
colData(sce)$mRNAsi <- mRNAsi
plotTSNE(sce, colour_by="mRNAsi") + fontsize

nogrids <- theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), axis.line = element_line(colour = "black"))
temp <- as.data.frame(reducedDim(sce, "TSNE"))

temp$col <- t(mRNAsi)
png("tSNE_mRNAsi.png", width = 1000, height = 1000)
ggplot(temp, aes(x=V1, y=V2)) + geom_point(size=3, aes(color=col)) + nogrids + scale_colour_gradient2(low = "blue", mid = "yellow", high = "red", midpoint = median(temp$col), space = "Lab", guide = "colourbar", aesthetics = "colour") + theme(legend.position="none") 
dev.off()

# violin plot, shown in Extended Data Figure 1n (right)
collist<-c("cluster_igraph_reordered", "mRNAsi")
temp <- colData(sce)[,collist]

temp$mRNAsi <- as.numeric(temp$mRNAsi)
temp <- as.data.frame(temp)

data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

png("violin_mRNAsi_reordered_Igraph.png", width = 1000, height = 1000)
ggplot(temp, aes(x=cluster_igraph_reordered, y=mRNAsi, fill=cluster_igraph_reordered)) + nogrids + geom_violin() + stat_summary(fun.data=data_summary) + theme(legend.position="none")
dev.off()

################################################################################

## ANOVA ANALYSES

## Statistical analysis for the differentiation scores in the individual clusters, for Fig.1 and Extended Data Fig.1
## First, make a data frame containing the cluster# and all of the differentiation scores for the individual cells.
## These scores are already calculated and saved as colData in the sce object. Thus, we can just retrieve the saved data:

## If you like to check a list of the colData in the sce object, use:
## colnames(colData(sce))

## Retrieve the cluster# and the differentiation scores. 
## as.data.frame is required otherwise the difscores object remains to be S4Vector, which is not compatible with ggplots2.
difscores <- colData(sce)[,c("cluster_igraph_reordered", "ecto.Z", "meso.Z", "endo.Z", "pluri.Z")]
difscores <- as.data.frame(difscores)
head(difscores)

## One example to draw a violin plot representing pluri.Z scores distribution in the individual clusters:
png("violin_pluriZ_cluster_igraph_reordered.png", width = 900, height = 1200)
ggplot(difscores, aes(x=cluster_igraph_reordered, y=pluri.Z, fill=cluster_igraph_reordered)) + geom_violin(scale = "width", trim = FALSE, lwd = 1.5) + theme_classic(base_size = 36) + labs(x="Clusters", y="PLURI score", title="PLURI score") + theme(plot.title = element_text(size=50, hjust = 0.5, face="bold"), axis.text.x = element_text(face="bold", size=32)) + theme(aspect.ratio = 1.33) + theme(legend.position="none") 
dev.off()

## One example to perform ANOVA:
pluri.model <- aov(pluri.Z~cluster_igraph_reordered, data=difscores)
summary(pluri.model) 
TukeyHSD(pluri.model, conf.level = 0.99)

################################################################################

## GENERATION OF MARKER LISTS 
## Find markers for the individual clusters
new.names <- rownames(sce)
dup.name <- new.names %in% new.names[duplicated(new.names)]
summary(dup.name)

new.names[dup.name] <- paste0(new.names, "_", rowData(sce)$ENSEMBL)[dup.name]
rownames(sce) <- new.names
temp <- duplicated(rownames(sce))
summary(temp)

markers <- findMarkers(sce, clusters=sce$cluster_igraph_reordered, assay.type="logcounts", get.spikes=FALSE, full.stats=FALSE)

for (i in 1:10){
  write.table(markers[[i]], file=paste(i, "BulkH9_markers_reordered_igraph_Cluster.tsv", sep = "_"), sep="\t", quote=FALSE, col.names=NA)
}
