################################################################################
# Orlando et al. Nature 2021
# R script used for endoderm markers in hepatocellular carcinoma vs healthy liver 
# written by: Dr. Juan Luis García-Rodríguez
# Date: May 10, 2020
################################################################################

# The Hepatocellular carcinoma data set is publicly available on NCBI, GEO accession: GSE103867
# The Human Healthy Liver data set is publicly available on NCBI, GEO accession: GSE115469  

################################################################################

# load libraries
library(DropletUtils)
library(Seurat)
library(dplyr)
library(ggplot2)
library(cowplot)
library(ggpubr)
library(tools)

################################################################################

## LOAD and ANALYZE THE COUNT MATRIX: Hepatocellular carcinoma

sce.data <- Read10X(data.dir = '~/folder_name/')
sce <- CreateSeuratObject(counts = sce.data, project = "ENDO_HCC", min.cells = 2, min.features = 200)
sce[["percent.mt"]] <- PercentageFeatureSet(sce, pattern = "^MT-")
VlnPlot(sce, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(sce, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
sce <- subset(sce, subset = nFeature_RNA > 200 & nFeature_RNA < 5500 & percent.mt <10 )
sce <- NormalizeData(sce, normalization.method = "LogNormalize", scale.factor = 10000)
sce <- FindVariableFeatures(sce, selection.method = "vst", nfeatures = 3000)
top10 <- head(VariableFeatures(sce), 10)
plot1 <- VariableFeaturePlot(sce)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
all.genes <- rownames(sce)
sce <- ScaleData(sce, vars.to.regress = "percent.mt")
sce <- RunPCA(sce, features = VariableFeatures(object = sce))
DimPlot(sce, reduction = "pca")
sce <- FindNeighbors(sce, dims = 1:30) 
sce <- FindClusters(sce, resolution = 0.5) 
head(Idents(sce), 5)
sce <- RunUMAP(sce, dims = 1:30)
DimPlot(sce, reduction = "umap")
sce <- RunTSNE(sce, dims = 1:30)
DimPlot(sce, reduction = "tsne")
sce.markers <- FindAllMarkers(sce, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
sce.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(sce.markers, "~/folder_name/sce_markers.csv")
saveRDS(sce, file = "~/folder_name/FILE.rds")
DefaultAssay(sce) <- "RNA"

################################################################################

## STUDY FEATURE LIST: Hepatocellular carcinoma 
# List file contains top 3 ENDO MARKERS: "RRP15", "CIRH1A", "HTR1D"

Common_Genelists <- file.path("~/List_Endo_top3")
umapTable <- as.data.frame(sce@reductions$tsne@cell.embeddings, stringsAsFactors = F)
plot_list <- list()
list_of_genlists <- list.files(Common_Genelists, pattern = ".txt")
plot_list[[1]] <- DimPlot(sce, reduction = "tsne", label = TRUE, order = TRUE) + NoLegend()
for (i in 1:length(list_of_genlists)) {
  cur_genelist_file <- file.path(Common_Genelists, list_of_genlists[i])
  genelist1 <- scan(cur_genelist_file, what = "", sep = "\t")
  genename <- file_path_sans_ext(list_of_genlists[i])
  sce <- AddModuleScore(object = sce, list(features = genelist1), ctrl = 25, name = genename)
  plot_list[[i + 1]] <- (FeaturePlot(object = sce, reduction = "tsne", features = paste(genename, "1", sep = ""), 
                                     cols = c("cyan", "red"), order = TRUE))
  dt_size <- dim(sce@meta.data)
  nCol <- dt_size[2]
  umapTable$tmp <- sce@meta.data[[nCol]]
  colnames(umapTable)[colnames(umapTable) == "tmp"] <- genename
}

################################################################################

## CALCULATE STEM CELL INDEX: Hepatocellular carcinoma

# Extract the counts as a matrix
temp <- sce@assays[["RNA"]]@counts
temp <- as.matrix(temp)
temp <- temp[!duplicated(row.names(temp)), ]
dim(temp)
write.table(temp, file = "~/folder_name/logcounts.txt", sep="\t", row.names=TRUE, col.names=TRUE)
# Add "gene_id" at the head of the "logcounts.txt". 

# Load the stem cell signature
w <- read.delim("~/folder_name/pcbc-stemsig.tsv", header=FALSE, row.names=1) %>% as.matrix() %>% drop()


## Reduces HUGO|POSITION gene IDs to just HUGO
f <- function( v ) unlist( lapply( strsplit( v, "\\|" ), "[[", 1 ) )


X <- read.delim("~/folder_name/logcounts.txt", as.is=TRUE, check.names=FALSE) %>%  ## Read the raw values
  filter( !grepl( "\\?", gene_id ) ) %>%      ## Drop genes with no mapping to HUGO
  mutate( gene_id = f( gene_id ) ) %>%        ## Clip gene ids to HUGO
  filter( gene_id %in% names(w) )         ## Reduce to the signature's gene set
## SLC35E2 has multiple entries with the same HUGO id
## Keep the first entry only
j <- grep( "SLC35E2", X[,1] )
if( length(j) > 1 )
  X <- X[-j[-1],]

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

write.table(cbind(s), file = "~/folder_name/StemScore.tsv", sep = "\t", quote = FALSE, col.names = FALSE)

mRNAsi <- read.delim("~/folder_name_name/StemScore.tsv", header=FALSE, row.names=1)
cellID <- colnames(sce)
mRNAsi <- mRNAsi$V2[match(cellID, rownames(mRNAsi))]
sce@meta.data$mRNAsi <- mRNAsi

# Selecting high vs low mRNAsi (Stem Cell Index score) at 0.5:
plot <- FeatureScatter(sce, "mRNAsi", "ENDO_Top3") + geom_rug(col="steelblue",alpha=0.1, size=1.5)
sce <- CellSelector(plot = plot, object = sce, ident = "High_mRNAsi")
sce <- CellSelector(plot = plot, object = sce, ident = "Low_mRNAsi")

Low_mRNAsi <- subset(sce, idents = c("Low_mRNAsi"))
Low_mRNAsi <- subset(sce, subset = mRNAsi < 0.5)
High_mRNAsi <- subset(sce, idents = c("High_mRNAsi"))
High_mRNAsi <- subset(sce, subset = mRNAsi > 0.5)
merged_High_Low_mRNAsi <- merge(High_mRNAsi, Low_mRNAsi) 

Low_mRNAsi <- subset(sce, idents = c("Low_mRNAsi"))
High_mRNAsi <- subset(sce, idents = c("High_mRNAsi"))
merged_High_Low_mRNAsi <- merge(High_mRNAsi, Low_mRNAsi) 

# Violin Plots shown in Extended Data Figure 5m:
ggpubr::annotate_figure(p = VlnPlot(merged_High_Low_mRNAsi, features = c("mRNAsi"), pt.size = 0, )) + geom_boxplot(width=0.1, color="black", alpha=0.2)  + rotate_x_text(angle = 45) + geom_hline(yintercept = mean(merged_High_Low_mRNAsi$mRNAsi), linetype = 2) + stat_compare_means(method = "t.test", label.y = 1600) + stat_compare_means(label = "p.format", method = "t.test", ref.group = "Low_mRNAsi")
ggsave("~/folder_name/ANOVA_mRNAsi_High_Low.pdf", width = 6, height = 4)

ggpubr::annotate_figure(p = VlnPlot(merged_High_Low_mRNAsi, features = c("ENDO_Top3"), pt.size = 0) + geom_boxplot(width=0.1, color="black", alpha=0.2)  + rotate_x_text(angle = 45) + geom_hline(yintercept = mean(merged_High_Low_mRNAsi$ENDO_Top31), linetype = 2) + stat_compare_means(method = "t.test", label.y = 1600) + stat_compare_means(label = "p.format", method = "t.test", ref.group = "Low_mRNAsi"))
ggsave("~/folder_name/ANOVA_ENDO_Top3_High_Low_mRNAsi.pdf", width = 6, height = 4)

################################################################################
################################################################################

## LOAD and ANALYZE THE COUNT MATRIX: Healthy Liver

countsData<-read.csv(file = "~/folder_name/GSE115469_Data.csv.gz", header = TRUE, sep = ",", row.names = 1)
Healthy_liver <- CreateSeuratObject(counts = countsData, project = "Healthy_liver", min.cells = 20, min.features = 200)
Healthy_liver[["percent.mt"]] <- PercentageFeatureSet(Healthy_liver, pattern = "^MT-")
VlnPlot(Healthy_liver, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)
plot1 <- FeatureScatter(Healthy_liver, feature1 = "nCount_RNA", feature2 = "percent.mt")
plot2 <- FeatureScatter(Healthy_liver, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
CombinePlots(plots = list(plot1, plot2))
Healthy_liver <- subset(Healthy_liver, subset = nFeature_RNA > 200 & nFeature_RNA < 6000 & percent.mt <10 )
Healthy_liver <- NormalizeData(Healthy_liver, normalization.method = "LogNormalize", scale.factor = 10000)
Healthy_liver <- FindVariableFeatures(Healthy_liver, selection.method = "vst", nfeatures = 3000)
top10 <- head(VariableFeatures(Healthy_liver), 10)
plot1 <- VariableFeaturePlot(Healthy_liver)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
CombinePlots(plots = list(plot1, plot2))
all.genes <- rownames(Healthy_liver)
Healthy_liver <- ScaleData(sce, vars.to.regress = "percent.mt")
Healthy_liver <- RunPCA(Healthy_liver, features = VariableFeatures(object = Healthy_liver))
DimPlot(Healthy_liver, reduction = "pca")
Healthy_liver <- FindNeighbors(Healthy_liver, dims = 1:30) 
Healthy_liver <- FindClusters(Healthy_liver, resolution = 0.5) 
head(Idents(Healthy_liver), 5)
Healthy_liver <- RunUMAP(Healthy_liver, dims = 1:30)
DimPlot(Healthy_liver, reduction = "umap")
Healthy_liver <- RunTSNE(Healthy_liver, dims = 1:30)
DimPlot(Healthy_liver, reduction = "tsne")
Healthy_liver.markers <- FindAllMarkers(Healthy_liver, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
Healthy_liver.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_logFC)
write.csv(Healthy_liver.markers, "~/folder_name/sce_markers.csv")
saveRDS(Healthy_liver, file = "~/folder_name/FILE.rds")
DefaultAssay(Healthy_liver) <- "RNA"

################################################################################

## STUDY FEATURE LISTS: Healthy Liver
# List file contains top 3 ENDO MARKERS: "RRP15", "CIRH1A", "HTR1D"

Common_Genelists <- file.path("~/List_Endo_top3")
umapTable <- as.data.frame(Healthy_liver@reductions$tsne@cell.embeddings, stringsAsFactors = F)
plot_list <- list()
list_of_genlists <- list.files(Common_Genelists, pattern = ".txt")
plot_list[[1]] <- DimPlot(Healthy_liver, reduction = "tsne", label = TRUE, order = TRUE) + NoLegend()
for (i in 1:length(list_of_genlists)) {
  cur_genelist_file <- file.path(Common_Genelists, list_of_genlists[i])
  genelist1 <- scan(cur_genelist_file, what = "", sep = "\t")
  genename <- file_path_sans_ext(list_of_genlists[i])
  Healthy_liver <- AddModuleScore(object = Healthy_liver, list(features = genelist1), ctrl = 25, name = genename)
  plot_list[[i + 1]] <- (FeaturePlot(object = Healthy_liver, reduction = "tsne", features = paste(genename, "1", sep = ""), 
                                     cols = c("cyan", "red"), order = TRUE))
  dt_size <- dim(Healthy_liver@meta.data)
  nCol <- dt_size[2]
  umapTable$tmp <- Healthy_liver@meta.data[[nCol]]
  colnames(umapTable)[colnames(umapTable) == "tmp"] <- genename
}

################################################################################

## CALCULATE STEM CELL INDEX: Healthy Liver

#Extract the counts and put them on a matrix.
temp <- Healthy_liver@assays[["RNA"]]@counts
temp <- as.matrix(temp)
temp <- temp[!duplicated(row.names(temp)), ]
dim(temp)
write.table(temp, file = "~/folder_name/logcounts.txt", sep="\t", row.names=TRUE, col.names=TRUE)
## Add "gene_id" at the head of the "logcounts.txt". 

#Load the stem cell signature
w <- read.delim("~/folder_name/pcbc-stemsig.tsv", header=FALSE, row.names=1) %>% as.matrix() %>% drop()

## Reduces HUGO|POSITION gene IDs to just HUGO
f <- function( v ) unlist( lapply( strsplit( v, "\\|" ), "[[", 1 ) )

X <- read.delim("~/folder_name/HCC_logcounts.txt", as.is=TRUE, check.names=FALSE) %>%  ## Read the raw values
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

write.table(cbind(s), file = "~/folder_name/StemScore.tsv", sep = "\t", quote = FALSE, col.names = FALSE)

mRNAsi <- read.delim("~/folder_name/StemScore.tsv", header=FALSE, row.names=1)
cellID <- colnames(Healthy_liver)
mRNAsi <- mRNAsi$V2[match(cellID, rownames(mRNAsi))]
Healthy_liver@meta.data$mRNAsi <- mRNAsi

# Selecting high vs low mRNAsi (Stem Cell Index score) at 0.5 score:
ENDO_Healthy_Liver <- Healthy_liver
plot <- FeatureScatter(ENDO_Healthy_Liver, "mRNAsi", "ENDO_Top31") + geom_rug(col="steelblue",alpha=0.1, size=1.5)
ENDO_Healthy_Liver <- CellSelector(plot = plot, object = ENDO_Healthy_Liver, ident = "High_mRNAsi")
ENDO_Healthy_Liver <- CellSelector(plot = plot, object = ENDO_Healthy_Liver, ident = "Low_mRNAsi")

Low_mRNAsi <- subset(ENDO_Healthy_Liver, idents = c("Low_mRNAsi"))
Low_mRNAsi <- subset(ENDO_Healthy_Liver, subset = mRNAsi < 0.5)
High_mRNAsi <- subset(ENDO_Healthy_Liver, idents = c("High_mRNAsi"))
High_mRNAsi <- subset(ENDO_Healthy_Liver, subset = mRNAsi > 0.5)
merged_High_Low_mRNAsi <- merge(High_mRNAsi, Low_mRNAsi) 

# Violin Plots shown in Extended Data Figure 5n:
ggpubr::annotate_figure(p = VlnPlot(merged_High_Low_mRNAsi, features = c("mRNAsi"), pt.size = 0, )) + geom_boxplot(width=0.1, color="black", alpha=0.2)  + rotate_x_text(angle = 45) + geom_hline(yintercept = mean(merged_High_Low_mRNAsi$mRNAsi), linetype = 2) + stat_compare_means(method = "t.test", label.y = 1600) + stat_compare_means(label = "p.format", method = "t.test", ref.group = "Low_mRNAsi")
ggsave("~/folder_name/ANOVA_mRNAsi_High_Low_mRNAsi.pdf", width = 6, height = 4)

ggpubr::annotate_figure(p = VlnPlot(merged_High_Low_mRNAsi, features = c("ENDO_Top3"), pt.size = 0) + geom_boxplot(width=0.1, color="black", alpha=0.2)  + rotate_x_text(angle = 45) + geom_hline(yintercept = mean(merged_High_Low_mRNAsi$ENDO_Top31), linetype = 2) + stat_compare_means(method = "t.test", label.y = 1600) + stat_compare_means(label = "p.format", method = "t.test", ref.group = "Low_mRNAsi"))
ggsave("~/folder_name/ANOVA_ENDO_Top3_High_Low_mRNAsi.pdf", width = 6, height = 4)
