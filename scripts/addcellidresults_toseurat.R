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

Seurat_file <- as.character(args[1])
doublet_file <- as.character(args[2])

print(Seurat_file)
print(doublet_file)

date = gsub("2022-","22",Sys.Date(),perl=TRUE);
date = gsub("-","",date);

seurat_object <- readRDS(Seurat_file)
doublet_object <- readRDS(doublet_file)
