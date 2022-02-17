.libPaths( c("/storage1/fs1/allegra.petti/Active/R_libs_scratch/RLibs_4.0.3",.libPaths()) )
library(SingleR)
library(SingleCellExperiment)
library(scuttle)
library(scran)
library(scater)
library(scRNAseq)
library(Seurat)
library(BiocParallel)
library(yaml)
library(Ckmeans.1d.dp)
library(stringr)
library(scSorter)
library(RColorBrewer)


args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 5) {
  args <- c("--help")
}

Seurat_rds <- as.character(args[1])
sample_name <- as.character(args[2])
scsorter_marker_list <- as.character(args[3])
marker_list_name <- as.character(args[4])


run_scsorter <- function(marker_list,seurat_obj,yaml_obj,def_wt=2){
  DefaultAssay(seurat_obj) <- 'RNA'
  annot_df_scsorter <- markerlist_to_df(marker_list=marker_list,seurat_obj=seurat_obj)
  expr_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000, verbose = F)
  topgenes <- head(VariableFeatures(expr_obj), 2000)
  expr = GetAssayData(expr_obj)
  topgene_filter = rowSums(as.matrix(expr)[topgenes, ]!=0) > ncol(expr)*.1
  topgenes = topgenes[topgene_filter]
  picked_genes = unique(c(annot_df_scsorter$Marker, topgenes))
  expr = expr[rownames(expr) %in% picked_genes, ]
  rts <- scSorter(expr, annot_df_scsorter,default_weight = def_wt)
  celltypes <- rts$Pred_Type
  names(celltypes) <- colnames(expr)
  #celltypes_mod <- paste0('nef.Malignant.',celltypes)
  #names(celltypes_mod) <- names(celltypes)
  return(celltypes)
}

scsorter_preds <- run_scsorter(marker_list = sel_scsorter_marker_list,seurat_obj = seurat_obj_malig,yaml_obj=yaml_obj)

saveRDS(scsorter_preds,paste0(sample_name,'_scsorter_preds_',marker_list_name,'.rds')
