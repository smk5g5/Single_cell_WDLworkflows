.libPaths( c("/storage1/fs1/allegra.petti/Active/R_libs_scratch/RLibs_4.0.3",.libPaths()) )
library(DoubletFinder)
library(scds)
library(scDblFinder)
library(scran)
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(Matrix)
library(DoubletCollection)
library(Seurat)
library(BiocParallel)
library(yaml)
library(Ckmeans.1d.dp)
library(stringr)
library(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 2) {
  args <- c("--help")
}

#This script will be updated to work with WDL

Seurat_file <- as.character(args[1])
doublet_file <- as.character(args[2])

print(Seurat_file)
print(doublet_file)

date = gsub("2022-","22",Sys.Date(),perl=TRUE);
date = gsub("-","",date);

seurat_object <- readRDS(Seurat_file)
doublet_object <- readRDS(doublet_file)


#seurat_object <- readRDS('/storage1/fs1/allegra.petti/Active/Users/khan.saad/WDL_pipelines/seurat_counts_to_anndata/0d1c2574-28a1-49e1-88ce-537ad2e57312/call-Get_seurat_counts/inputs/-1835829567/Cycling.SCT.PCA.UMAP.TSNE.CLUST.211216.rds')