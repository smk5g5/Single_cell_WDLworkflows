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
if(length(args) < 4) {
  args <- c("--help")
}

Seurat_rds <- as.character(args[1])
sample_name <- as.character(args[2])
scsorter_marker_list_rds <- as.character(args[3])
marker_list_name <- as.character(args[4])

markerlist_to_df <- function(marker_list,seurat_obj){
  sub_marker_list <- list()
  df <- data.frame(Type=character(),Marker=character())
  for(i in names(marker_list)){
    sub_marker_list[[i]] <- intersect(marker_list[[i]],rownames(seurat_obj))
  }
  for(i in names(sub_marker_list)){
    celltype <- rep(i,length(sub_marker_list[[i]]))
    de <- list(Type=celltype,Marker=sub_marker_list[[i]])
    df = rbind(df,de,stringsAsFactors=FALSE)
  }
  return(df)
}

run_scsorter <- function(marker_list,seurat_obj,def_wt=2){
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

seurat_obj <- readRDS(Seurat_rds)

scsorter_marker_list <- readRDS(scsorter_marker_list_rds)

sel_scsorter_marker_list <- lapply(scsorter_marker_list,intersect,x=rownames(seurat_obj))

scsorter_preds <- run_scsorter(marker_list = sel_scsorter_marker_list,seurat_obj = seurat_obj)

saveRDS(scsorter_preds,paste0(sample_name,'_scsorter_preds_',marker_list_name,'.rds'))
