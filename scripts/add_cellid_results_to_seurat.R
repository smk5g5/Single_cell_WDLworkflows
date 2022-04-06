.libPaths( c("/storage1/fs1/allegra.petti/Active/R_libs_scratch/RLibs_4.0.3",.libPaths()) )
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 4) {
  args <- c("--help")
}

date = gsub("2022-","22",Sys.Date(),perl=TRUE);
date = gsub("-","",date);


Seurat_file <- as.character(args[1])
cellid_pred_file <- as.character(args[2])
cellid_mat_file <- as.character(args[3])
marker_file_reference_name <- as.character(args[4])

seurat_object <- readRDS(Seurat_file)
cellid_pred <- readRDS(cellid_pred_file)
cellid_mat <- readRDS(cellid_mat_file)


Add_cellid_scores_to_seuratassay <- function(cellid_mat,seurat_object,refname){
	cellid_assay <- CreateAssayObject(data = t(as.matrix(cellid_mat)))
    seurat_object[[sprintf("%s.CELLID", refname)]] <- cellid_assay
  }
return(seurat_object)
}

Add_cellid_preds_to_seuratmeta <- function(cellid_pred,seurat_object,refname){
  index <- match(Cells(seurat_object),names(cellid_pred))
  seurat_object[[sprintf("CELLID_results_%s",refname)]] <- cellid_pred[index]
  return(seurat_object)  
}

seurat_object <- Add_cellid_scores_to_seuratassay(cellid_mat=cellid_mat,seurat_object=seurat_object,refname=marker_file_reference_name)

seurat_object <- Add_cellid_preds_to_seuratmeta(cellid_pred=cellid_pred,seurat_object=seurat_object,refname=marker_file_reference_name)


saveRDS(seurat_object, file = paste0('Cycling.SCT.PCA.UMAP.TSNE.CLUST.',".Cellid_results.",date,".rds"))
