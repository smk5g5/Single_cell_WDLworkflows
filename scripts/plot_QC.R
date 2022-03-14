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

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 4) {
  args <- c("--help")
}

seurat_rds <- as.character(args[1])
organism <- as.character(args[2])
project_name <- as.character(args[3])
pre_post <- as.character(args[4])


date = gsub("2022-","22",Sys.Date(),perl=TRUE);
date = gsub("-","",date);


scrna_GEX <- readRDS(seurat_rds)

plot_qc_metrics <- function(seurat_obj,filename,mt_pat) {
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


if(organism=='human'){
plot_qc_metrics(seurat_obj=scrna_GEX,filename=sprintf("%s.%s_filtering_QC.%s.pdf",project_name,pre_post,date),mt_pat="^MT-")
} else{
plot_qc_metrics(seurat_obj=scrna_GEX,filename=sprintf("%s.%s_filtering_QC.%s.pdf",project_name ,pre_post, date),mt_pat="^mt-")
}