################################################################################
# Orlando et al. Nature 2021
# R script used for scRNAseq transformed V1H9 (trans, t-PSC) versus healthy H9 hESC data analysis 
# written by: Dr. Mio Nakanishi
# Date: October 3, 2019
################################################################################

# The H9 count matrix file is publicly available on NCBI, GEO accession: GSE126022 
# Originally published: Nakanishi et al. Cell 2019. PMID: 30982595

# The V1H9 count matrix file is publicly available on NCBI, GEO accession: GSE133669

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
library(akima)

################################################################################

## H9V1: LOAD COUNT MATRIX FROM GEO AND CREATE STEM CELL EXPERIMENT

# set working directory
# download transformed V1H9 count matrix from on NCBI, GEO accession: GSE133669

# load data
V1H9_GEO <- read.table("v1H9_counts.txt", header=TRUE,sep="\t")
# remove any potentially duplicated rows
V1H9_GEO <- V1H9_GEO[!duplicated(V1H9_GEO$X), ]
rownames(H9_GEO) <- H9_GEO[,1]
V1H9_GEO[,1] <- NULL # remove column X
V1H9_GEO[,1] <- NULL # remove column X.1

# set-up stem cell experiment (sce)
all.counts <- V1H9_GEO
all.counts <- as.matrix(all.counts)
cell_labels <- colnames(V1H9_GEO)
genes <- rownames(V1H9_GEO)

sce <- SingleCellExperiment(list(counts=all.counts),
                            colData=DataFrame(cell_labels),
                            rowData=DataFrame(genes),
                            metadata=list(study="V1H9_GEO"))

# QC FOR V1H9
dim(sce)

is.mito <- grepl("^MT-", rowData(sce)$Symbol)
summary(is.mito)

sce <- calculateQCMetrics(sce, feature_controls=list(Mt=is.mito)) 

libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE)
feature.drop <- isOutlier(sce$total_features_by_counts, nmads=3, type="lower", log=TRUE)
mito.drop <- isOutlier(sce$pct_counts_Mt, nmads=3, type="higher", log=TRUE)
sce <- sce[,!(libsize.drop | feature.drop | mito.drop)]
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), ByMito=sum(mito.drop), Remaining=ncol(sce))

sce.v1H9 <- sce

ave.counts.v1H9 <- calcAverage(sce.v1H9, use_size_factors=FALSE)
rowData(sce.v1H9)$ave.count <- ave.counts.v1H9
to.keep.v1H9 <- ave.counts.v1H9 > 0
sce.v1H9 <- sce.v1H9[to.keep.v1H9,]
summary(to.keep.v1H9)

clusters.v1H9 <- quickCluster(sce.v1H9, min.mean=0.1, method="igraph")
sce.v1H9 <- computeSumFactors(sce.v1H9, cluster=clusters.v1H9, min.mean=0.1)
summary(sizeFactors(sce.v1H9))

sce.v1H9 <- normalize(sce.v1H9)

var.fit.nospike.v1H9 <- trendVar(sce.v1H9, parametric=TRUE, use.spikes=FALSE, loess.args=list(span=0.2))
var.dec.nospike.v1H9 <- decomposeVar(sce.v1H9, var.fit.nospike.v1H9)

var.dec.nospike.v1H9$Symbol <- rowData(sce.v1H9)$Symbol
var.dec.nospike.v1H9 <- var.dec.nospike.v1H9[order(var.dec.nospike.v1H9$bio, decreasing=TRUE),]
head(var.dec.nospike.v1H9)

################################################################################

## LOAD COUNT MATRIX FROM GEO AND CREATE STEM CELL EXPERIMENT

# download H9 count matrix from on NCBI, GEO accession: GSE126022
# load data
H9_GEO <- read.table("GSE126022_bulk_NCADpos_H9_counts.txt", header=TRUE,sep="\t")
# Cells were removed if they were present in the NCAD+ sorted population (refer to Figure 1 script)

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
dim(sce)

is.mito <- grepl("^MT-", rowData(sce)$Symbol)
summary(is.mito)

sce <- calculateQCMetrics(sce, feature_controls=list(Mt=is.mito)) 

libsize.drop <- isOutlier(sce$total_counts, nmads=3, type="lower", log=TRUE)
feature.drop <- isOutlier(sce$total_features_by_counts, nmads=3, type="lower", log=TRUE)
mito.drop <- isOutlier(sce$pct_counts_Mt, nmads=3, type="higher", log=TRUE)
sce <- sce[,!(libsize.drop | feature.drop | mito.drop)]
data.frame(ByLibSize=sum(libsize.drop), ByFeature=sum(feature.drop), ByMito=sum(mito.drop), Remaining=ncol(sce))

BulkSex <- read.delim(file="sex_annotation_Bulk_H1_H9.txt", header=TRUE)
BulkSex$Cell_Barcode <- paste(sep="",BulkSex$Cell_Barcode,"-1")
SexAnnotation<- BulkSex$Sex[match(colnames(sce), BulkSex$Cell_Barcode)]
sce$Sex <- SexAnnotation
sce.H9 <- filter(sce, sce$Sex == "female")
dim(sce.H9)

ave.counts.H9 <- calcAverage(sce.H9, use_size_factors=FALSE)
rowData(sce.H9)$ave.count <- ave.counts.H9
to.keep.H9 <- ave.counts.H9 > 0
sce.H9 <- sce.H9[to.keep.H9,]
summary(to.keep.H9)

clusters.H9 <- quickCluster(sce.H9, min.mean=0.1, method="igraph")
sce.H9 <- computeSumFactors(sce.H9, cluster=clusters.H9, min.mean=0.1)
summary(sizeFactors(sce.H9))

sce.H9 <- normalize(sce.H9)

colnames(sce.H9) <- gsub(pattern="1", replacement="2", colnames(sce.H9))

var.fit.nospike.H9 <- trendVar(sce.H9, parametric=TRUE, use.spikes=FALSE, loess.args=list(span=0.2))
var.dec.nospike.H9 <- decomposeVar(sce.H9, var.fit.nospike.H9)

var.dec.nospike.H9$Symbol <- rowData(sce.H9)$Symbol
var.dec.nospike.H9 <- var.dec.nospike.H9[order(var.dec.nospike.H9$bio, decreasing=TRUE),]
head(var.dec.nospike.H9)

################################################################################

## MNN BATCH CORRECTION

universe <- intersect(rownames(var.dec.nospike.v1H9), rownames(var.dec.nospike.H9))
mean.bio <- (var.dec.nospike.v1H9[universe,"bio"] + var.dec.nospike.H9[universe,"bio"])/2
chosen <- universe[mean.bio > 0]
length(chosen)

rescaled <- multiBatchNorm(sce.v1H9[universe,], sce.H9[universe,])
rescaled.v1H9 <- rescaled[[1]]
rescaled.H9 <- rescaled[[2]]

set.seed(100) 
original <- list(
  v1H9=logcounts(rescaled.v1H9)[chosen,],
  H9=logcounts(rescaled.H9)[chosen,]
)

B.v1H9<-as.matrix(logcounts(rescaled.v1H9)[chosen,])
B.H9<-as.matrix(logcounts(rescaled.H9)[chosen,])

mnn.out3 <- mnnCorrect(B.v1H9, B.H9, k=20, cos.norm.out=FALSE)

omat <- do.call(cbind, original)
sce <- SingleCellExperiment(list(logcounts=omat))
reducedDim(sce, "MNN") <- mnn.out3$corrected
sce$Batch <- as.character(mnn.out3$batch)
sce

dim(mnn.out3$corrected[[1]])
dim(mnn.out3$corrected[[2]])

mnn.exp <-cbind(mnn.out3$corrected[[1]],mnn.out3$corrected[[2]])
grch38 <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")

################################################################################

## FIGURE 2A: PLURI score on trans V1H9 versus healthy H9

temp.markers <- scan("List_Gifford_PluriMarkers.txt", what="", sep="\t")
temp.markers <- getBM(attributes ="ensembl_gene_id", filters = "external_gene_name", values = temp.markers, mart = grch38, uniqueRows=F)
temp.markers <- intersect(rownames(mnn.exp), temp.markers$ensembl_gene_id)
temp.vals <- as.matrix((mnn.exp)[temp.markers,,drop=FALSE])
temp.sig <- colMeans(temp.vals)/nrow(temp.vals) - colMeans(mnn.exp)/nrow(mnn.exp)
temp.sig <- as.matrix(temp.sig)
temp.z <- temp.sig/sd(temp.sig) - mean(temp.sig)/sd(temp.sig)
temp <- as.data.frame(temp.z)
temp$Batch <- as.factor(c(rep("v1H9", 2427),rep("H9",2186)))

png("violin_Gifford_Pluri_v1H9_H9.png", width = 900, height = 1200)
ggplot(temp, aes(x=Batch, y=V1, fill=Batch)) + geom_violin(scale = "width", trim = FALSE, lwd = 1.5) + theme_classic(base_size = 36) + labs(x="Clusters", y="ECTO scores", title="ECTO score") + theme(plot.title = element_text(size=58, hjust = 0.5, face="bold"), axis.text.x = element_text(face="bold", size=32)) + theme(aspect.ratio = 1.33) + theme(legend.position="none") +geom_signif(comparisons=list(c("v1H9","H9")), map_signif_level=TRUE, y_position=6, textsize=12)
dev.off()

################################################################################

## FIGURE 2B: Stem Cell Index score on trans V1H9 versus healthy H9

temp <- mnn.exp[!duplicated(row.names(mnn.exp)), ]
rownames(temp) <- rowData(sce.H9)$Symbol[match(rownames(temp),rownames(sce.H9))]
temp <- temp[!duplicated(row.names(temp)), ]
write.table(temp, file = "v1H9_H9_mnnCorrect.txt", sep="\t", row.names=TRUE, col.names=TRUE)

w <- read.delim("pcbc-stemsig.tsv", header=FALSE, row.names=1) %>% as.matrix() %>% drop()

## Reduces HUGO|POSITION gene IDs to just HUGO
f <- function( v ) unlist( lapply( strsplit( v, "\\|" ), "[[", 1 ) )

X <- read.delim("v1H9_H9_mnnCorrect.txt", as.is=TRUE, check.names=FALSE, sep="\t") %>%  ## Read the raw values
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

write.table(cbind(s), file = "v1H9_H9_scRNA-seq_mnnCorrect_StemScore.tsv", sep = "\t", quote = FALSE, col.names = FALSE)

mRNAsi <- read.delim("v1H9_H9_scRNA-seq_mnnCorrect_StemScore.tsv", header=FALSE, row.names=1)
cellID <- colnames(mnn.exp)
mRNAsi <- mRNAsi$V2[match(cellID, rownames(mRNAsi))]
mRNAsi <- t(mRNAsi)

temp <- as.data.frame(t(mRNAsi))
temp$Batch <- as.factor(c(rep("v1H9", 2427),rep("H9",2186)))

png("violin_mRNAsi_v1H9_H9.png", width = 900, height = 1200)
ggplot(temp, aes(x=Batch, y=V1, fill=Batch)) + geom_violin(scale = "width", trim = FALSE, lwd = 1.5) + theme_classic(base_size = 36) + labs(x="Clusters", y="Stem Cell Index", title="Stem Cell Index") + theme(plot.title = element_text(size=58, hjust = 0.5, face="bold"), axis.text.x = element_text(face="bold", size=32)) + theme(aspect.ratio = 1.33) + theme(legend.position="none") +geom_signif(comparisons=list(c("v1H9","H9")), map_signif_level=TRUE, y_position=1.2, textsize=12)
dev.off()

################################################################################

## tSNE with courier, trans V1H9 and healthy H9 hESCs
## Shown in Figure 2c

V1.markers <- scan("List_V1H9_UPmarkers.txt", what="", sep="\t")
H9.markers <- scan("List_V1H9_DOWNmarkers.txt", what="", sep="\t")

V1.markers <- intersect(rownames(sce),  V1.markers)
V1.vals <- as.matrix(logcounts(sce)[V1.markers,,drop=FALSE])
all.vals <- as.matrix(logcounts(sce))
V1.sig <- colMeans(V1.vals)/nrow(V1.vals) - colMeans(all.vals)/nrow(all.vals)
V1.sig <- as.matrix(V1.sig)
V1.z <- V1.sig/sd(V1.sig) - mean(V1.sig)/sd(V1.sig)
colData(sce)$V1.Z <- V1.z

H9.markers <- intersect(rownames(sce), H9.markers)
H9.vals <- as.matrix(logcounts(sce)[H9.markers,,drop=FALSE])
all.vals <- as.matrix(logcounts(sce))
H9.sig <- colMeans(H9.vals)/nrow(H9.vals) - colMeans(all.vals)/nrow(all.vals)
H9.sig <- as.matrix(H9.sig)
H9.z <- H9.sig/sd(H9.sig) - mean(H9.sig)/sd(H9.sig)
colData(sce)$H9.Z <- H9.z

temp <- as.data.frame(reducedDim(sce, "TSNE"))
temp$V1 <- V1.z
V1_H9.z <- V1.z - H9.z
colData(sce)$V1_H9.Z <- V1_H9.z
temp$V1_H9 <- V1_H9.z

temp.grid <- interp(as.data.frame(temp)$V1, as.data.frame(temp)$V2, as.data.frame(temp)$V1)
temp.grid2 <- expand.grid(x=temp.grid$x, y=temp.grid$y)
temp.grid2$z <- as.vector(temp.grid$z)

png("BulkH9_tSNE_V1_SSEA3pos.adjPVal0.005.up_grid.png", width = 1000, height = 1000)
ggplot(as.data.frame(temp), aes(x=V1, y=V2))  + geom_tile(data=na.omit(temp.grid2), aes(x=x, y=y, fill=z)) + scale_fill_gradient2(low = "white", mid = "white", high = "red") + geom_point(size=2, color="gray25")+ theme_classic(base_size = 40) + labs(x="tSNE 1", y="tSNE 2", title="v1H9 signature")  + theme(plot.title = element_text(size=50, hjust = 0.5, face="bold"))+ theme(aspect.ratio = 1) 
dev.off()

################################################################################

## FIGURE 2C: PLURI score on trans V1H9 versus healthy H9

# Meso Score
temp.markers <- scan("List_Gifford_MesoMarkers.txt", what="", sep="\t")
temp.markers <- getBM(attributes ="ensembl_gene_id", filters = "external_gene_name", values = temp.markers, mart = grch38, uniqueRows=F)
temp.markers <- intersect(rownames(mnn.exp), temp.markers$ensembl_gene_id)
temp.vals <- as.matrix((mnn.exp)[temp.markers,,drop=FALSE])
temp.sig <- colMeans(temp.vals)/nrow(temp.vals) - colMeans(mnn.exp)/nrow(mnn.exp)
temp.sig <- as.matrix(temp.sig)
temp.z <- temp.sig/sd(temp.sig) - mean(temp.sig)/sd(temp.sig)
temp <- as.data.frame(temp.z)
temp$Batch <- as.factor(c(rep("v1H9", 2427),rep("H9",2186)))

png("violin_Gifford_Meso_v1H9_H9.png", width = 900, height = 1200)
ggplot(temp, aes(x=Batch, y=V1, fill=Batch)) + geom_violin(scale = "width", trim = FALSE, lwd = 1.5) + theme_classic(base_size = 36) + labs(x="Clusters", y="ECTO scores", title="ECTO score") + theme(plot.title = element_text(size=58, hjust = 0.5, face="bold"), axis.text.x = element_text(face="bold", size=32)) + theme(aspect.ratio = 1.33) + theme(legend.position="none") +geom_signif(comparisons=list(c("v1H9","H9")), map_signif_level=TRUE, y_position=6, textsize=12)
dev.off()

# Ecto Score
temp.markers <- scan("List_Gifford_EctoMarkers.txt", what="", sep="\t")
temp.markers <- getBM(attributes ="ensembl_gene_id", filters = "external_gene_name", values = temp.markers, mart = grch38, uniqueRows=F)
temp.markers <- intersect(rownames(mnn.exp), temp.markers$ensembl_gene_id)
temp.vals <- as.matrix((mnn.exp)[temp.markers,,drop=FALSE])
temp.sig <- colMeans(temp.vals)/nrow(temp.vals) - colMeans(mnn.exp)/nrow(mnn.exp)
temp.sig <- as.matrix(temp.sig)
temp.z <- temp.sig/sd(temp.sig) - mean(temp.sig)/sd(temp.sig)
temp <- as.data.frame(temp.z)
temp$Batch <- as.factor(c(rep("v1H9", 2427),rep("H9",2186)))

png("violin_Gifford_Ecto_v1H9_H9.png", width = 900, height = 1200)
ggplot(temp, aes(x=Batch, y=V1, fill=Batch)) + geom_violin(scale = "width", trim = FALSE, lwd = 1.5) + theme_classic(base_size = 36) + labs(x="Clusters", y="ECTO scores", title="ECTO score") + theme(plot.title = element_text(size=58, hjust = 0.5, face="bold"), axis.text.x = element_text(face="bold", size=32)) + theme(aspect.ratio = 1.33) + theme(legend.position="none") +geom_signif(comparisons=list(c("v1H9","H9")), map_signif_level=TRUE, y_position=6, textsize=12)
dev.off()

# Endo Score
temp.markers <- scan("List_Gifford_EndoMarkers.txt", what="", sep="\t")
temp.markers <- getBM(attributes ="ensembl_gene_id", filters = "external_gene_name", values = temp.markers, mart = grch38, uniqueRows=F)
temp.markers <- intersect(rownames(mnn.exp), temp.markers$ensembl_gene_id)
temp.vals <- as.matrix((mnn.exp)[temp.markers,,drop=FALSE])
temp.sig <- colMeans(temp.vals)/nrow(temp.vals) - colMeans(mnn.exp)/nrow(mnn.exp)
temp.sig <- as.matrix(temp.sig)
temp.z <- temp.sig/sd(temp.sig) - mean(temp.sig)/sd(temp.sig)
temp <- as.data.frame(temp.z)
temp$Batch <- as.factor(c(rep("v1H9", 2427),rep("H9",2186)))

png("violin_Gifford_Endo_v1H9_H9.png", width = 900, height = 1200)
ggplot(temp, aes(x=Batch, y=V1, fill=Batch)) + geom_violin(scale = "width", trim = FALSE, lwd = 1.5) + theme_classic(base_size = 36) + labs(x="Clusters", y="ECTO scores", title="ECTO score") + theme(plot.title = element_text(size=58, hjust = 0.5, face="bold"), axis.text.x = element_text(face="bold", size=32)) + theme(aspect.ratio = 1.33) + theme(legend.position="none") +geom_signif(comparisons=list(c("v1H9","H9")), map_signif_level=TRUE, y_position=6, textsize=12)
dev.off()
