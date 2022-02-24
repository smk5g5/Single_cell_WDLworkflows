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
# if(length(args) < 4) {
#   args <- c("--help")
# }
args_len <- length(args)
Seurat_rds <- args[1]
sample_name <- args[2]
type <- args[3]


date = gsub("2022-","22",Sys.Date(),perl=TRUE);
date = gsub("-","",date);


Seurat_obj <- readRDS(Seurat_rds)

for(i in 4:args_len){
	rds_name <- args[i]
	ref_name <- gsub('.rds','',unlist(str_split(basename(rds_name),pattern='_scsorter_preds_',simplify=T))[2])
	ref_results <- readRDS(rds_name)
	index <- match(Cells(Seurat_obj),names(ref_results))
	Seurat_obj[[paste0('scSorter.',type,'.',ref_name)]] <- ref_results[index]
}

saveRDS(Seurat_obj,paste0(sample_name,'.',type,'.seurat_scsorter_mergedpreds.',date,'.rds'))