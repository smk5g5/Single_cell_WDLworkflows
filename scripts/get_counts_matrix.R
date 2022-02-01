.libPaths( c("/storage1/fs1/allegra.petti/Active/R_libs_scratch/RLibs_4.0.3",.libPaths()) )
library(Seurat)
library(dplyr)


args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1) {
  args <- c("--help")
}

Seurat_file <- as.character(args[1])

seurat_object <- readRDS(Seurat_file)

