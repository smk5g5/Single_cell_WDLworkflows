#!/usr/local/bin/Rscript

.libPaths( c("/storage1/fs1/allegra.petti/Active/R_libs_scratch/RLibs_4.0.3",.libPaths()) )
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(Matrix)
#source('./scripts/Plot_QC_scrnaseq.R')

plot_qc_metrics <- function(seurat_obj,filename,mt_pat="^MT-") {
seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
seurat_obj$mitoRatio <- PercentageFeatureSet(object = seurat_obj, pattern = mt_pat)
seurat_obj$mitoRatio <- seurat_obj@meta.data$mitoRatio / 100
metadata <- seurat_obj@meta.data
metadata$cells <- rownames(metadata)

metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

Ncells <- metadata %>% 
  ggplot(aes(x=Sample, fill=Sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

cell_dens <- metadata %>% 
  ggplot(aes(color=Sample, x=nUMI, fill= Sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)+ggtitle("Cell density")

genes_per_cell <- metadata %>% 
  ggplot(aes(color=Sample, x=nGene, fill= Sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300) + ggtitle("Genes per cell")

box_genespercell <- metadata %>% 
  ggplot(aes(x=Sample, y=log10(nGene), fill=Sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

umibygenesmito <- metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~Sample) + ggtitle("UMI vs Genes colored by mito ratio")


mitoratio_dens <- metadata %>% 
  ggplot(aes(color=Sample, x=mitoRatio, fill=Sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2) + ggtitle("mitoratio density")

genes_per_umi <- metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = Sample, fill=Sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8) + ggtitle("Genes per UMI")

pdf(filename)
print(Ncells)
print(cell_dens)
print(genes_per_cell)
print(box_genespercell)
print(umibygenesmito)
print(mitoratio_dens)
print(genes_per_umi)
dev.off()
}

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

# ##################################################################################
# # This bit is only for testing purposes for the rscript within wdl would be disabled
# # in the main workflow
set.seed(100)
random_sample_of_cells = sample(Cells(scrna_GEX),2500)
#select 10% of all cells randomly for testing the script.
scrna_GEX <- subset(scrna_GEX,cells=random_sample_of_cells)
###################################################################################


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

png(sprintf("%s.VlnPlot.Filtered.nCount.25Kmax.control.%s.png", sample_name, date), width = 10, height = 10, units="in", res=300)
vln <- VlnPlot(object = scrna_GEX, features = "nCount_RNA", pt.size=0, group.by="Sample", y.max=25000)
print(vln)
dev.off();
#
png(sprintf("%s.VlnPlot.Filtered.nFeature.control.%s.png", sample_name,date), width = 10, height = 10, units="in", res=300)
vln <- VlnPlot(object = scrna_GEX, features = "nFeature_RNA", pt.size=0, group.by="Sample")
print(vln)
dev.off()

saveRDS(scrna_GEX, file = sprintf("%s.MergedFilteredSeuratObject.%s.rds", sample_name, date))