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
# if(length(args) < 4) {
#   args <- c("--help")
# }

Seurat_file <- as.character(args[1])
singleR_tsv_file <- as.character(args[2])
singleR_files <- as.character(args[3:length(args)])

date = gsub("2022-","22",Sys.Date(),perl=TRUE);
date = gsub("-","",date);



input_df <- read.table(singleR_tsv_file,sep="\t",header=FALSE)
colnames(input_df) <- c('Reference_name','label_column_name','reference_rds')


single_R_preds <- list()

for(i in 1:nrow(input_df)){

	
}