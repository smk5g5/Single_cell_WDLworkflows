#!/usr/local/bin/Rscript

.libPaths( c("/storage1/fs1/allegra.petti/Active/R_libs_scratch/RLibs_4.0.3",.libPaths()) )
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(Matrix)
source('./scripts/Plot_QC_scrnaseq.R')


args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2) {
  args <- c("--help")
}

Seurat_10x_directory <- as.character(args[1])
sample_name <- as.character(args[2])

genome <- "GRCh38";
date = gsub("2021-","21",Sys.Date(),perl=TRUE);
date = gsub("-","",date);

data.10x  <- Read10X(data.dir = Seurat_10x_directory);
scrna_GEX = CreateSeuratObject(counts = data.10x, min.cells=10, min.features=100, project=sample_name);
scrna_GEX$Sample = sample_name

mito.genes <- grep(pattern = "^MT-", x = rownames(x = scrna_GEX), value = TRUE,ignore.case = TRUE);
percent.mito <- Matrix::colSums(x = GetAssayData(object = scrna_GEX, slot = 'counts')[mito.genes, ]) / Matrix::colSums(x = GetAssayData(object = scrna_GEX, slot = 'counts'));
scrna_GEX[['percent.mito']] <- percent.mito;

ribo.genes <- grep(pattern = "^RP[SL][[:digit:]]", x = rownames(x = scrna_GEX), value = TRUE,ignore.case = TRUE);
percent.ribo <- Matrix::colSums(x = GetAssayData(object = scrna_GEX, slot = 'counts')[ribo.genes, ]) / Matrix::colSums(x = GetAssayData(object = scrna_GEX, slot = 'counts'));
scrna_GEX[['percent.ribo']] <- percent.ribo;

print("Making violin plots...");
#
png(sprintf("%s.VlnPlot.%s.png", sample_name, date), width = 13, height = 6, units="in", res=300);
vln <- VlnPlot(object = scrna_GEX, features = c("percent.mito", "percent.ribo"), pt.size=0, ncol = 2, group.by="Sample");
print(vln);
dev.off();

png(sprintf("%s.VlnPlot.nCount.25Kmax.%s.png", sample_name,date), width = 10, height = 10, units="in", res=300)
vln <- VlnPlot(object = scrna_GEX, features = "nCount_RNA", pt.size=0, group.by="Sample", y.max=25000)
print(vln)
dev.off();

png(sprintf("%s.VlnPlot.nFeature.%s.png", sample_name,date), width = 10, height = 10, units="in", res=300)
vln <- VlnPlot(object = scrna_GEX, features = "nFeature_RNA", pt.size=0, group.by="Sample")
print(vln)
dev.off()

print("Making scatter plots...")
pdf(sprintf("%s.Scatter1.%s.pdf", sample_name, date), width = 8, height = 6, useDingbats=FALSE);
scatter <- FeatureScatter(object = scrna_GEX, feature1 = "nCount_RNA", feature2 = "percent.mito", pt.size=0.1,group.by = 'Sample')
print(scatter);
dev.off();
pdf(sprintf("%s.Scatter2.%s.pdf", sample_name, date), width = 8, height = 6, useDingbats=FALSE);
scatter <- FeatureScatter(object = scrna_GEX, feature1 = "nCount_RNA", feature2 = "percent.ribo", pt.size=0.1,group.by = 'Sample')
print(scatter);
dev.off();
pdf(sprintf("%s.Scatter3.%s.pdf", sample_name, date), width = 8, height = 6, useDingbats=FALSE);
scatter <- FeatureScatter(object = scrna_GEX, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.1,group.by = 'Sample')
print(scatter);
dev.off();

saveRDS(scrna_GEX, file = sprintf("%s.unfilteredSeuratObject_GEX.%s.rds",sample_name , date))

plot_qc_metrics(seurat_obj=scrna_GEX,filename=sprintf("%s.Pre_filtering_QC.%s.pdf",sample_name , date))

mc.hi = 0.05;
m <- median(scrna_GEX@meta.data$nFeature_RNA)
s <- sd(scrna_GEX@meta.data$nFeature_RNA)
Feature95 <- quantile(scrna_GEX@meta.data$nFeature_RNA, 0.95);
Feature05 <- quantile(scrna_GEX@meta.data$nFeature_RNA, 0.05);
m1 <- mean(scrna_GEX@meta.data$nCount_RNA)
s1 <- sd(scrna_GEX@meta.data$nCount_RNA)
Count93 <- quantile(scrna_GEX@meta.data$nCount_RNA, 0.93);
scrna_GEX <- subset(x = scrna_GEX, subset = nFeature_RNA > 700 & nCount_RNA < Count93 & percent.mito < mc.hi)

plot_qc_metrics(seurat_obj=scrna_GEX,filename=sprintf("%s.Post_filtering_QC.%s.pdf",sample_name , date))

png(sprintf("%s.VlnPlot.Filtered.nCount.25Kmax.control.png", sample_name, date), width = 10, height = 10, units="in", res=300)
vln <- VlnPlot(object = scrna_GEX, features = "nCount_RNA", pt.size=0, group.by="Sample", y.max=25000)
print(vln)
dev.off();
#
png(sprintf("%s.VlnPlot.Filtered.nFeature.%s.control.png", sample_name, , date), width = 10, height = 10, units="in", res=300)
vln <- VlnPlot(object = scrna_GEX, features = "nFeature_RNA", pt.size=0, group.by="Sample")
print(vln)
dev.off()

saveRDS(scrna_GEX, file = sprintf("%s.MergedFilteredSeuratObject.%s.rds", sample_name, date))