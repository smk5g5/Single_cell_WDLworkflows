library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(Matrix)
library(RColorBrewer)
library(ggthemes)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 6) {
  args <- c("--help")
}

seurat_loc <- as.character(args[1])
sub_col <- as.character(args[2])
ident_names <- as.character(args[3])
inverse <- as.character(args[4])
# output_suffix <- as.character(args[5])
sample_name <- as.character(args[5])
# output.stats <- as.character(args[4])
# output_meta <- as.character(args[5])

seurat_obj <- readRDS(seurat_loc)

print(names(seurat_obj@meta.data))

date = gsub("2023-","23",Sys.Date(),perl=TRUE);
date = gsub("-","",date);

DefaultAssay(seurat_obj) <- "RNA"
# print('Does the error happen inside subset_renormalize_recluster?')
Idents(seurat_obj) <- sub_col
if(inverse=='TRUE'){
  scrna_GEX <- subset(seurat_obj,idents=ident_names,invert = TRUE)
}else{
  scrna_GEX <- subset(seurat_obj,idents=ident_names)
}

saveRDS(scrna_GEX,paste0('Cycling.SCT.PCA.UMAP.TSNE.CLUST.',sample_name,'.doublet_removed.230511.rds'))
