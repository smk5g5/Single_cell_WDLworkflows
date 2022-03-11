#!/usr/local/bin/Rscript

.libPaths( c("/storage1/fs1/allegra.petti/Active/R_libs_scratch/RLibs_4.0.3",.libPaths()) )
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(Matrix)
library("sctransform");
library("dplyr");
library("RColorBrewer");
library("ggthemes");
library("ggplot2");
library("cowplot");
library(tidyverse)

library(ggplot2)
library(gridExtra)
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
if(length(args) < 4) {
  args <- c("--help")
}

input_tsv_file <- as.character(args[1])
organism <- as.character(args[2])
project_name <- as.character(args[3])
gene_lists <- as.character(args[4])

if(organism=='human'){
genome <- "GRCh38";
}
else{
  genome <- "GRCm38";
}
date = gsub("2022-","22",Sys.Date(),perl=TRUE);
date = gsub("-","",date);

input_df <- read.table(input_tsv_file,sep="\t",header=FALSE)
colnames(input_df) <- c('Sample','cellranger_10x_directory')

# #input_tsv_file looks like this for example
# KO1 /storage1/fs1/allegra.petti/Active/GBM/Stegh/SAMPLES/KO1/filtered_feature_bc_matrix
# KO2 /storage1/fs1/allegra.petti/Active/GBM/Stegh/SAMPLES/KO2/filtered_feature_bc_matrix

cluster.res=0.7; #default cluster resolution
nPC=50; #default number of principal components
#final number of prinicipal component will be determined by the out
#put of jackstraw analysis.

matrix.dirs = list();

for(i in 1:nrow(input_df)){
  matrix.dirs[i] = input_df$cellranger_10x_directory[i]
}

control = "Cycling" #default might need to be changed in the next iteration
matrix.dirs <- unlist(matrix.dirs[]);
names(matrix.dirs) <- input_df$Sample

sample_names <- input_df$Sample

data.10x = list();
scrna_GEX.list = list();

for (i in 1:length(matrix.dirs)) {
 print(i);
 data.10x[[i]] <- Read10X(data.dir = matrix.dirs[i]);
 scrna_GEX.list[[i]] = CreateSeuratObject(counts = data.10x[[i]], min.cells=10, min.features=100, project=project_name);
 scrna_GEX.list[[i]][["Sample"]] = sample_names[i];
}

scrna_GEX <- merge(x=scrna_GEX.list[[1]],y=c(scrna_GEX.list[[2:length(sample_names)]]),add.cell.ids = sample_names,project=project_name)

mito.genes <- grep(pattern = "^MT-", x = rownames(x = scrna_GEX), value = TRUE,ignore.case = TRUE);
percent.mito <- Matrix::colSums(x = GetAssayData(object = scrna_GEX, slot = 'counts')[mito.genes, ]) / Matrix::colSums(x = GetAssayData(object = scrna_GEX, slot = 'counts'));
scrna_GEX[['percent.mito']] <- percent.mito;

ribo.genes <- grep(pattern = "^RP[SL][[:digit:]]", x = rownames(x = scrna_GEX), value = TRUE,ignore.case = TRUE);
percent.ribo <- Matrix::colSums(x = GetAssayData(object = scrna_GEX, slot = 'counts')[ribo.genes, ]) / Matrix::colSums(x = GetAssayData(object = scrna_GEX, slot = 'counts'));
scrna_GEX[['percent.ribo']] <- percent.ribo;

print("Making violin plots...");
#
png(sprintf("%s.VlnPlot.%s.png", project_name, date), width = 13, height = 6, units="in", res=300);
vln <- VlnPlot(object = scrna_GEX, features = c("percent.mito", "percent.ribo"), pt.size=0, ncol = 2, group.by="Sample");
print(vln);
dev.off();

png(sprintf("%s.VlnPlot.nCount.25Kmax.%s.png", project_name,date), width = 10, height = 10, units="in", res=300)
vln <- VlnPlot(object = scrna_GEX, features = "nCount_RNA", pt.size=0, group.by="Sample", y.max=25000)
print(vln)
dev.off();

png(sprintf("%s.VlnPlot.nFeature.%s.png", project_name,date), width = 10, height = 10, units="in", res=300)
vln <- VlnPlot(object = scrna_GEX, features = "nFeature_RNA", pt.size=0, group.by="Sample")
print(vln)
dev.off()

print("Making scatter plots...")
pdf(sprintf("%s.Scatter1.%s.pdf", project_name, date), width = 8, height = 6, useDingbats=FALSE);
scatter <- FeatureScatter(object = scrna_GEX, feature1 = "nCount_RNA", feature2 = "percent.mito", pt.size=0.1,group.by = 'Sample')
print(scatter);
dev.off();
pdf(sprintf("%s.Scatter2.%s.pdf", project_name, date), width = 8, height = 6, useDingbats=FALSE);
scatter <- FeatureScatter(object = scrna_GEX, feature1 = "nCount_RNA", feature2 = "percent.ribo", pt.size=0.1,group.by = 'Sample')
print(scatter);
dev.off();
pdf(sprintf("%s.Scatter3.%s.pdf", project_name, date), width = 8, height = 6, useDingbats=FALSE);
scatter <- FeatureScatter(object = scrna_GEX, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", pt.size=0.1,group.by = 'Sample')
print(scatter);
dev.off();

saveRDS(scrna_GEX, file = sprintf("MergedSeuratObject_GEX.%s.rds",date))

plot_qc_metrics(seurat_obj=scrna_GEX,filename=sprintf("%s.Pre_filtering_QC.%s.pdf",project_name , date))

mc.hi = 0.05;

if(organism=='human'){
mc.hi = 0.1;
}
else{
mc.hi = 0.05;
}

for (i in 1:length(x = scrna_GEX.list)) {
     m <- median(scrna_GEX.list[[i]]@meta.data$nFeature_RNA)
     s <- sd(scrna_GEX.list[[i]]@meta.data$nFeature_RNA)
     Feature95 <- quantile(scrna_GEX.list[[i]]@meta.data$nFeature_RNA, 0.95);
     Feature05 <- quantile(scrna_GEX.list[[i]]@meta.data$nFeature_RNA, 0.05);
     mito.genes <- grep(pattern = "^MT-", x = rownames(x = scrna_GEX.list[[i]]), value = TRUE,ignore.case = TRUE);
     percent.mito <- Matrix::colSums(x = GetAssayData(object = scrna_GEX.list[[i]], slot = 'counts')[mito.genes, ]) / Matrix::colSums(x = GetAssayData(object = scrna_GEX.list[[i]], slot = 'counts'));
     scrna_GEX.list[[i]][['percent.mito']] <- percent.mito;
     m1 <- mean(scrna_GEX.list[[i]]@meta.data$nCount_RNA)
     s1 <- sd(scrna_GEX.list[[i]]@meta.data$nCount_RNA)
     Count93 <- quantile(scrna_GEX.list[[i]]@meta.data$nCount_RNA, 0.93);
     scrna_GEX.list[[i]] <- subset(x = scrna_GEX.list[[i]], subset = nFeature_RNA > 700 & nCount_RNA < Count93 & percent.mito < mc.hi)
}

scrna_GEX <- merge(x=scrna_GEX.list[[1]],y=c(scrna_GEX.list[[2:length(sample_names)]]),add.cell.ids = sample_names,project=project_name)

plot_qc_metrics(seurat_obj=scrna_GEX,filename=sprintf("%s.Post_filtering_QC.%s.pdf",project_name, date))

png(sprintf("%s.VlnPlot.Filtered.nCount.25Kmax.control.%s.png", project_name, date), width = 10, height = 10, units="in", res=300)
vln <- VlnPlot(object = scrna_GEX, features = "nCount_RNA", pt.size=0, group.by="Sample", y.max=25000)
print(vln)
dev.off();
#
png(sprintf("%s.VlnPlot.Filtered.nFeature.control.%s.png", project_name,date), width = 10, height = 10, units="in", res=300)
vln <- VlnPlot(object = scrna_GEX, features = "nFeature_RNA", pt.size=0, group.by="Sample")
print(vln)
dev.off()

saveRDS(scrna_GEX, file = sprintf("%s.MergedFilteredSeuratObject.%s.rds", project_name, date))

png(sprintf("%s/VlnPlot.Filtered.nCount.25Kmax.%s.%s.png", output.stats, control, date), width = 10, height = 10, units="in", res=300)
vln <- VlnPlot(object = scrna_GEX, features = "nCount_RNA", pt.size=0, group.by="Batch", y.max=25000)
print(vln)
dev.off();

png(sprintf("%s/VlnPlot.Filtered.nFeature.%s.%s.png", output.stats, control, date), width = 10, height = 10, units="in", res=300)
vln <- VlnPlot(object = scrna_GEX, features = "nFeature_RNA", pt.size=0, group.by="Batch")
print(vln)
dev.off()


if(organism=='human'){
cell.cycle.tirosh <- read.table("/storage1/fs1/allegra.petti/Active/10xGenomics/key.gene.lists/CellCycleTirosh.txt", sep='\t', header=TRUE);
s.genes = cell.cycle.tirosh$`Gene.Symbol`[which(cell.cycle.tirosh$List == "G1/S")];
g2m.genes = cell.cycle.tirosh$`Gene.Symbol`[which(cell.cycle.tirosh$List == "G2/M")];
}
else{
cell.cycle.tirosh <- read.table("/storage1/fs1/allegra.petti/Active/10xGenomics/key.gene.lists/CellCycleTirosh_mouse.txt", sep='\t', header=FALSE);
s.genes = cell.cycle.tirosh$V2[which(cell.cycle.tirosh$V1 == "G1/S")];
g2m.genes = cell.cycle.tirosh$V2[which(cell.cycle.tirosh$V1 == "G2/M")];
}

scrna_GEX <- CellCycleScoring(object=scrna_GEX, s.features=s.genes, g2m.features=g2m.genes, set.ident=FALSE)

scrna_GEX <- NormalizeData(object = scrna_GEX, normalization.method = "LogNormalize", scale.factor = 1e6); # 1e6 is new as of 1/8/20

scrna_GEX <- FindVariableFeatures(object = scrna_GEX, selection.method = 'vst', mean.cutoff = c(0.1,8), dispersion.cutoff = c(1, Inf))
print(paste("Number of Variable Features: ",length(x = VariableFeatures(object = scrna_GEX))));
#
VG.file = sprintf("%s.variableGenes.%s.pdf",control, date);
pdf(VG.file, useDingbats=FALSE)
vg <- VariableFeaturePlot(scrna_GEX)
print(vg);
dev.off()

if (control == "Cycling") { # This removes all signal associated with the cell cycle
 scrna_GEX <- ScaleData(object = scrna_GEX, features = rownames(x = scrna_GEX), vars.to.regress = c("S.Score","G2M.Score"), display.progress=FALSE);
} else if (control == "CyclingRB") {
 scrna_GEX <- ScaleData(object = scrna_GEX, features = rownames(x = scrna_GEX), vars.to.regress = c("S.Score","G2M.Score","percent.ribo"), display.progress=FALSE);
} else if (control == "CyclingDiff") {
 scrna_GEX <- ScaleData(object = scrna_GEX, features = rownames(x = scrna_GEX), vars.to.regress = c("CC.Difference"), display.progress=FALSE);
} else {
 scrna_GEX <- ScaleData(object = scrna_GEX, features = rownames(x = scrna_GEX), display.progress=FALSE);
}

saveRDS(scrna_GEX, file = sprintf("VST.%s.rds", date))

scrna_GEX <- RunPCA(object = scrna_GEX, npcs = 50, verbose = FALSE);
elbow <- ElbowPlot(object = scrna_GEX)
jpeg(sprintf("PCA%d.elbow.%s.%s.jpg", nPC, control, date), width = 6, height = 8, units="in", res=300);
print(elbow);
dev.off();
scrna_GEX <- RunUMAP(object = scrna_GEX, reduction = "pca", dims = 1:50)
scrna_GEX <- RunTSNE(object = scrna_GEX, reduction = "pca", dims = 1:50)

scrna_GEX <- JackStraw(object = scrna_GEX, num.replicate = 100, dims=nPC)
scrna_GEX <- ScoreJackStraw(object = scrna_GEX, dims = 1:nPC)
jpeg(sprintf("PCA.jackstraw.%s.%s.jpg", output.stats, control, date), width = 10, height = 6, units="in", res=300);
js <- JackStrawPlot(object = scrna_GEX, dims = 1:nPC)
print(js);
dev.off();

pc.pval <- scrna_GEX@reductions$pca@jackstraw@overall.p.values
print(pc.pval);
write.table(pc.pval, file=sprintf("PCA.jackstraw.scores.%s.%s.xls", output.stats, control, date), quote=FALSE, sep='\t', col.names=TRUE);

nPC=length( pc.pval[,'Score'][pc.pval[,'Score'] <= 0.05]) #redefine nPCs based on number of significant prinicipal components in jackstraw plot

scrna_GEX <- RunPCA(object = scrna_GEX, npcs = nPC, verbose = FALSE);
elbow <- ElbowPlot(object = scrna_GEX)
jpeg(sprintf("PCA%d.elbow.%s.%s.jpg", nPC, control, date), width = 6, height = 8, units="in", res=300);
print(elbow);
dev.off();
scrna_GEX <- RunUMAP(object = scrna_GEX, reduction = "pca", dims = 1:nPC)
scrna_GEX <- RunTSNE(object = scrna_GEX, reduction = "pca", dims = 1:nPC)

scrna_GEX <- JackStraw(object = scrna_GEX, num.replicate = 100, dims=nPC)
scrna_GEX <- ScoreJackStraw(object = scrna_GEX, dims = 1:nPC)
jpeg(sprintf("PCA.jackstraw.%s.%s.jpg", output.stats, control, date), width = 10, height = 6, units="in", res=300);
js <- JackStrawPlot(object = scrna_GEX, dims = 1:nPC)
print(js);
dev.off();


jpeg(sprintf("UMAP.%s.%s.jpg", output.stats, control, date), width = 10, height = 8, units="in", res=300);
p2 <- DimPlot(object = scrna_GEX, reduction = "umap", group.by = "Sample", pt.size=0.1)
print(p2);
dev.off();

print ("VizDimLoadings Running...");
jpeg(sprintf("VizDimLoadings.%s.%s.jpg", output.stats, control, date), width = 8, height = 30, units="in", res=300);
vdl <- VizDimLoadings(object = scrna_GEX, dims = 1:3)
print(vdl);
dev.off();

print ("ProjectDim Running...");
scrna_GEX <- ProjectDim(object = scrna_GEX)



print ("DimHeatmap Running...");
jpeg(sprintf("PCA.heatmap.top.%s.%s.jpg", output.stats, control, date), width = 8.5, height = 11, units="in", res=300);
hm <- DimHeatmap(object = scrna_GEX, dims = 1, cells = 500, balanced = TRUE);
print(hm);
dev.off();

jpeg(sprintf("PCA.heatmap.multi.%s.%s.jpg", output.stats, control, date), width = 8.5, height = 24, units="in", res=300);
hm.multi <- DimHeatmap(object = scrna_GEX, dims = 1:10, cells = 500, balanced = TRUE);
print(hm.multi);
dev.off();


saveRDS(scrna_GEX, file = sprintf("%s.SCT.PCA.UMAP.TSNE.%s.rds",control, date))

scrna_GEX <- FindNeighbors(object=scrna_GEX, dims=1:nPC);
scrna_GEX <- FindClusters(object=scrna_GEX, resolution=cluster.res);
scrna_GEX[[sprintf("ClusterNames_%.1f_%dPC", cluster.res, nPC)]] <- Idents(object = scrna_GEX)

DEGs <- FindAllMarkers(object=scrna_GEX); # output is a matrix!
write.table(DEGs, file=sprintf("DEGs.Wilcox.PCA.%d.%s.%s.xls", nPC, control, date), quote=FALSE, sep="\t", row.names=FALSE) # must save cluster-specific marker genes
saveRDS(scrna_GEX, file = sprintf("%s/%s.SCT.PCA.UMAP.TSNE.CLUST.%s.rds", control, date))

control <- 'Cycling'
n.graph = length(unique(scrna_GEX[[sprintf("ClusterNames_%.1f_%dPC",cluster.res, nPC)]][,1]));
rainbow.colors = rainbow(n.graph, s=0.6, v=0.9);

jpeg(sprintf("UMAP.clusters.%d.%.1f.%s.%s.jpg", nPC, cluster.res, control, date), width = 10, height = 8, units="in", res=300);
p2 <- DimPlot(object = scrna_GEX, reduction = "umap", group.by = sprintf("ClusterNames_%.1f_%dPC",cluster.res, nPC), cols = rainbow.colors, pt.size=0.1) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
print(p2);
dev.off();

jpeg(sprintf("UMAP.clusters.labeled.%d.%.1f.%s.%s.jpg", nPC, cluster.res, control, date), width = 10, height = 8, units="in", res=300);
p2 <- DimPlot(object = scrna_GEX, reduction = "umap", group.by = sprintf("ClusterNames_%.1f_%dPC",cluster.res, nPC), cols = rainbow.colors, pt.size=0.1, label=TRUE) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
print(p2);
dev.off();

feature.pal = rev(colorRampPalette(brewer.pal(11,"Spectral"))(50));

mito.genes <- grep(pattern = "^MT-", x = rownames(x = scrna_GEX), value = TRUE,ignore.case = TRUE);
percent.mito <- Matrix::colSums(x = GetAssayData(object = scrna_GEX, slot = 'counts')[mito.genes, ]) / Matrix::colSums(x = GetAssayData(object = scrna_GEX, slot = 'counts'));
scrna_GEX[['percent.mito']] <- percent.mito;

# ribosomal genes
ribo.genes <- grep(pattern = "^RP[SL][[:digit:]]", x = rownames(x = scrna_GEX), value = TRUE,ignore.case = TRUE);
percent.ribo <- Matrix::colSums(x = GetAssayData(object = scrna_GEX, slot = 'counts')[ribo.genes, ]) / Matrix::colSums(x = GetAssayData(object = scrna_GEX, slot = 'counts'));
scrna_GEX[['percent.ribo']] <- percent.ribo;

scrna_GEX[['percent.ribo']] <- percent.ribo;

print("Making additional UMAP plots");
# color UMAP plots by parameters of interest:
  print ("color by UMI");
  jpeg(sprintf("umap.%d.%.1f.colorby.UMI.%s.%s.jpg", nPC, cluster.res, control, date), width = 10, height = 8,  units="in", res=300);
  fp2 <- FeaturePlot(object = scrna_GEX, features = c("nCount_RNA"), cols = feature.pal, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  print(fp2);
  dev.off();

print ("color by % mito");
  jpeg(sprintf("umap.%d.%.1f.colorby.MC.%s.%s.jpg", nPC, cluster.res, control, date), width = 10, height = 8, units="in", res=300);
  fp2 <- FeaturePlot(object = scrna_GEX, features = c("percent.mito"), cols = feature.pal, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  print(fp2);
  dev.off();

print ("color by % RB");
  jpeg(sprintf("umap.%d.%.1f.colorby.RB.%s.%s.jpg", nPC, cluster.res, control, date), width = 10, height = 8, units="in", res=300);
  fp2 <- FeaturePlot(object = scrna_GEX, features = c("percent.ribo"), cols = feature.pal, pt.size=0.1, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  print(fp2);
  dev.off();

print ("color by nGene");
  jpeg(sprintf("umap.%d.%.1f.colorby.nGene.%s.%s.jpg", nPC, cluster.res, control, date), width = 10, height = 8, units="in", res=300);
  fp2 <- FeaturePlot(object = scrna_GEX, features = c("nFeature_RNA"), cols = feature.pal, pt.size=0.1, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  print(fp2);
 dev.off();

print ("color by Phase");
  jpeg(sprintf("umap.%d.%.1f.colorby.Phase.%s.%s.jpg", nPC, cluster.res, control, date), width = 10, height = 8, units="in", res=300);
  phase.colors = ptol_pal()(3)
  umapplot <- DimPlot(object = scrna_GEX, cols=phase.colors, group.by="Phase", pt.size=0.1, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank())
  print(umapplot);
  dev.off();
 
print ("color UMAP by Principal Components");
  jpeg(sprintf("UMAP.%d.%.1f.colorby.PCs.%s.%s.jpg", nPC, cluster.res, control, date), width = 12, height = 6, units="in", res=100);
  redblue=c("blue","gray","red");
  fp1 <- FeaturePlot(object = scrna_GEX, features = 'PC_1', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  fp2 <- FeaturePlot(object = scrna_GEX, features = 'PC_2', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  fp3 <- FeaturePlot(object = scrna_GEX, features = 'PC_3', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  fp4 <- FeaturePlot(object = scrna_GEX, features = 'PC_4', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  fp5 <- FeaturePlot(object = scrna_GEX, features = 'PC_5', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  fp6 <- FeaturePlot(object = scrna_GEX, features = 'PC_6', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  fp7 <- FeaturePlot(object = scrna_GEX, features = 'PC_7', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  fp8 <- FeaturePlot(object = scrna_GEX, features = 'PC_8', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  fp9 <- FeaturePlot(object = scrna_GEX, features = 'PC_9', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  fp10 <- FeaturePlot(object = scrna_GEX, features = 'PC_10', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
 print(plot_grid(fp1, fp2, fp3, fp4, fp5, fp6, fp7, fp8, fp9, fp10));
 # print(plot_grid(fp1, fp2, fp3, fp4, fp5, fp6, fp7, fp8, fp9));
  dev.off();

saveRDS(scrna_GEX, file = sprintf("%s.SCT.PCA.UMAP.TSNE.CLUST.%s.rds",control,date))

#must have column header containing List,Name 
# where List is the broad celltype and Name contains the marker gene for that celltype
#sep="," separator has to be comma for this

glists.raw <- read.table(gene_lists, sep=",",row.names=NULL,header=TRUE,as.is=TRUE); # gene lists
glists <- glists.raw[which(glists.raw$Name %in% rownames(scrna_GEX)),]; # filtered genelists

for(i in unique(glists$List)) {
  j=gsub(" ","_",i);
  j=gsub("/","_",j);
  sub_glists <- subset(glists,List==i)
  genesToPlot = sub_glists$Name
  ng = length(genesToPlot); # number of genes
  outfile = sprintf("umap.%s.%s.%s.pdf",j, control, date);
  if(ng==1){
  fp <- FeaturePlot(object = scrna_GEX, features = genesToPlot, cols = c("gray","red"), ncol=2, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  ggsave(outfile,plot=fp,width = 5, height = 5)
  }
  else{
    gplotlist=list();
    for(k in genesToPlot){
      gplotlist[[k]] <- FeaturePlot(object = scrna_GEX, features = genesToPlot, cols = c("gray","red"), ncol=2, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank()) + ggtitle(sprintf("%s (%s)",k,j))
    }
    ml <- marrangeGrob(gplotlist, nrow=2, ncol=2)
    ggsave(outfile,ml,width=10, height=10)
  }
}
