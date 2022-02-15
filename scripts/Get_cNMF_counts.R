.libPaths( c("/storage1/fs1/allegra.petti/Active/R_libs_scratch/RLibs_4.0.3",.libPaths()) )
library(Seurat)
library(BiocParallel)
library(yaml)
library(Ckmeans.1d.dp)
library(stringr)
library("RColorBrewer")
library("ggthemes");
library("ggplot2");
library("cowplot");
library(MASS)


date = gsub("2021-","21",Sys.Date(),perl=TRUE);
date = gsub("-","",date);

args = commandArgs(trailingOnly=TRUE);
print (length(args))

if(length(args) != 2)
{
  stop("\n too few or too many arguments provided \n")
}


rds_input <- args[1]
transposed_count_out <- args[2]

seurat_object <- readRDS(rds_input)
# seurat_object <- readRDS('/storage1/fs1/alberthkim/Active/users/khan.saad/Scenic_iteration2/schwann_newlabels.rds')
mycounts_mat <- seurat_object@assays$RNA@counts


#remove genes where all cells have zero count
#This is a cNMF requirements that
#cells which have zero count for a all genes should be removed
#as well as genes which have zero count for all cells should be removed

mycounts_mat <- t(as.matrix(mycounts_mat))



rownames_to_rem <- rownames(mycounts_mat[rowSums(mycounts_mat)==0,])

colnames_to_rem  <- colnames(mycounts_mat[,colSums(mycounts_mat)==0])

if( (length(rownames_to_rem)>0)&(length(colnames_to_rem)>0) ){
  	sel_rownames <- setdiff(rownames(mycounts_mat),rownames_to_rem)
  	sel_colnames <- setdiff(colnames(mycounts_mat),colnames_to_rem)
  	sel_matrix <-  mycounts_mat[sel_rownames,sel_colnames]
  	write.table(sel_matrix,sprintf("%s_RNA_assay_counts_matrix.%s.txt",prefix,date),sep = "\t",row.names=T,col.names=NA)
  }else if( (length(rownames_to_rem)>0)&(length(colnames_to_rem)==0) ){
	sel_rownames <- setdiff(rownames(mycounts_mat),rownames_to_rem)
  	sel_matrix <-  mycounts_mat[sel_rownames,]
	write.table(sel_matrix,sprintf("%s_RNA_assay_counts_matrix.%s.txt",prefix,date),sep = "\t",row.names=T,col.names=NA)
 }else if( (length(rownames_to_rem)==0)&(length(colnames_to_rem)>0) ){
  	sel_colnames <- setdiff(colnames(mycounts_mat),colnames_to_rem)
  	sel_matrix <-  mycounts_mat[,sel_colnames]
  	write.table(sel_matrix,sprintf("%s_RNA_assay_counts_matrix.%s.txt",prefix,date),sep = "\t",row.names=T,col.names=NA)
 }else{
 	write.table(mycounts_mat,sprintf("%s_RNA_assay_counts_matrix.%s.txt",prefix,date),sep = "\t",row.names=T,col.names=NA)
 }




# sel_rownames <- setdiff(rownames(mycounts_mat),rownames_to_rem)

# sel_matrix <-  mycounts_mat[sel_rownames,]


# write.table(mycounts_mat,file = transposed_count_out,sep = "\t",col.names = TRUE,row.names = TRUE)

# prepare_cmd = 'python ../cnmf.py prepare --output-dir %s --name %s -c %s -k %s --n-iter %d --total-workers 1 --seed %d --numgenes %d --beta-loss frobenius' % (output_directory, run_name, countfn, K, numiter, seed, numhvgenes)

# print(head(mycounts_mat))
# print(head(as.matrix(mycounts_mat)))
# write.table(mycounts_mat,file = transposed_count_out,sep = "\t",col.names = TRUE,row.names = TRUE)
# #fwrite(as.matrix(mycounts_mat),file = transposed_count_out,sep = "\t",row.names=T,col.names=T)
