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
singleR_reference <- as.character(args[2])
reference_name =  as.numeric(args[3])
label_column_name =  as.numeric(args[4])


seurat_obj <- readRDS(Seurat_rds)
singleR_obj <- readRDS(singleR_reference)

singleR_preds <- SingleR(test = as.SingleCellExperiment(seurat_obj),
                             ref = singleR_obj,
                             labels = singleR_obj[[label_column_name]],
                             BPPARAM=MulticoreParam())

saveRDS(object=singleR_preds,file=sprintf("%s_singleR_preds.rds", reference_name))
