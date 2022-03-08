.libPaths( c("/storage1/fs1/allegra.petti/Active/R_libs_scratch/RLibs_4.0.3",.libPaths()) )
library(Seurat)
library(dplyr)
library(MASS)



args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1) {
  args <- c("--help")
}

Seurat_file <- as.character(args[1])
prefix <- as.character(args[2])

date = gsub("2021-","21",Sys.Date(),perl=TRUE);
date = gsub("-","",date);

seurat_object <- readRDS(Seurat_file)

DefaultAssay(seurat_object) <- 'RNA'

counts_matrix <- as.matrix(seurat_obj@assays$RNA@counts)

counts_matrix <- t(counts_matrix)

write.matrix(counts_matrix,sprintf("%s_RNA_assay_counts_matrix_transposed.%s.txt",prefix,date),sep = "\t")


write.table(sub_metadata_df, file="sub_metadata_df_scrna.txt", row.names=T, col.names=T,na='',sep="\t")

write_counts_matrix <- function(seurat_obj,prefix,date) {

DefaultAssay(seurat_object) <- 'RNA'

counts_matrix <- as.matrix(seurat_object@assays$RNA@counts)

counts_matrix <- t(counts_matrix)

write.table(counts_matrix,sprintf("%s_RNA_assay_counts_matrix_transposed.%s.txt",'schwann',date),sep = "\t",row.names=T,col.names=NA)

}




# seurat_object <- macs

# DefaultAssay(seurat_object) <- 'RNA'

# write.table(counts_matrix,sprintf("%s_RNA_assay_counts_matrix_transposed.%s.txt",'macs_scrna',date),sep = "\t")

# write.table(counts_matrix, sprintf("%s_RNA_assay_counts_matrix_transposed.%s.txt",'macs_scrna',date),sep = "\t",row.names=TRUE,col.names=TRUE)


# "orig.ident"	"cluster_functionalannot"

# metadata_df <- seurat_object@meta.data

# sub_metadata_df <- metadata_df[c('orig.ident','newlabel')]

# colnames(sub_metadata_df) <- c("orig.ident","cluster_functionalannot")

# write.table(sub_metadata_df, file=sprintf("%s_RNA_assay_counts_matrix_transposed.%s.txt",'macs_scrna',date), row.names=T, col.names=NA,sep="\t")


# get_counts_matrix <- function(seurat_object,prefix,date){

# DefaultAssay(seurat_object) <- 'RNA'

# counts_matrix <- as.matrix(seurat_object@assays$RNA@counts)

# counts_matrix <- t(counts_matrix)

# write.table(counts_matrix,sprintf("%s_RNA_assay_counts_matrix_transposed.%s.txt",prefix,date),sep = "\t",row.names=T,col.names=NA)

# metadata_df <- seurat_object@meta.data

# sub_metadata_df <- metadata_df[c('orig.ident','newlabel')]

# colnames(sub_metadata_df) <- c("orig.ident","cluster_functionalannot")

# write.table(sub_metadata_df, file=sprintf("%s_sub_metadata_df_scrna.%s.txt",prefix,date), row.names=T, col.names=NA,sep="\t")
# }



# get_counts_matrix(seurat_object,prefix,date)

