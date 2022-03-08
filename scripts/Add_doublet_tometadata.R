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

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 3) {
  args <- c("--help")
}

Seurat_file <- as.character(args[1])
doublet_file <- as.character(args[2])
outfile_seurat <- as.character(args[3])

seurat_object <- readRDS(Seurat_file)
doublet_object <- readRDS(doublet_file)

add_doublet_predictions_to_seurat <- function(seurat_object,doublet_object){
  doublet_object[['doublet_results']] <- FindDoublets(score.list=doublet_object$doublet_scores,rate=0.08)
  for(i in names(doublet_object[['doublet_results']])){
    sel_index <- doublet_object[['doublet_results']][[i]]
    doublet_cells <- doublet_object$cellnames[sel_index]
    non_doublet_cells <- setdiff(doublet_object$cellnames,doublet_cells)
    doublet_cells_pred <- rep('yes',length(doublet_cells))
    names(doublet_cells_pred) <- doublet_cells
    non_doublet_cells_pred <- rep('no',length(non_doublet_cells))
    names(non_doublet_cells_pred) <- non_doublet_cells
    doublet_preds <- c(doublet_cells_pred,non_doublet_cells_pred)
    sel_cells <- intersect(Cells(seurat_object),names(doublet_preds))
    sel_doublet_preds <- doublet_preds[sel_cells]
    index <- match(Cells(seurat_object),names(sel_doublet_preds))
    seurat_object[[sprintf("doublet_results_%s",i)]] <- sel_doublet_preds[index]
  }
  for(i in names(doublet_object[['doublet_scores']])){
    doubscores <- doublet_object[['doublet_scores']][[i]]
    names(doubscores) <- doublet_object$cellnames
    sel_cells <- intersect(Cells(seurat_object),names(doubscores))
    sel_doubscores <- doubscores[sel_cells]
    index <- match(Cells(seurat_object),names(sel_doubscores))
    seurat_object[[sprintf("doublet_scores_%s",i)]] <- sel_doubscores[index]
  }
  return(seurat_object)
}


seurat_object <- add_doublet_predictions_to_seurat(seurat_object=seurat_object,doublet_object=doublet_object)
saveRDS(object = seurat_object,file=outfile_seurat)
