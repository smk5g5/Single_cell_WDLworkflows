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
# if(length(args) < 4) {
#   args <- c("--help")
# }



Seurat_file <- as.character(args[1])
input_tsv_file <- as.character(args[2])
project_name <- as.character(args[3])
split_by <- as.character(args[4])
doublet_files <- as.character(args[5:length(args)])

date = gsub("2022-","22",Sys.Date(),perl=TRUE);
date = gsub("-","",date);


add_doublet_predictions_to_seurat_singlesample_mod <- function(seurat_object,doublet_object,seurat_cellname_prefix){
  # doublet_object[['doublet_results']] <- FindDoublets(score.list=doublet_object$doublet_scores,rate=0.08)
  mydoublet_pred_df <- list()
  mydoublet_score_df <- list()
  for(i in names(doublet_object[['doublet_results']])){
    sel_index <- doublet_object[['doublet_results']][[i]]
    print(head(sel_index))
    new_cellnames <- paste0(seurat_cellname_prefix,'_',seurat_cellname_prefix,'_',doublet_object$cellnames)
    #This might be variable but based on how we pass sample names in 
    #seurat_multisample_input_pca_clustering.R script we get something like this as
    #barcode name with sample name repeated in the following format
    #this might not work unless seurat_multisample_input_pca_clustering.R is used for
    # generating the initial merged seurat object
    #use at your own risk!
    #rightnow the cell barcode name will be like
    #seurat_cellname_prefix_seurat_cellname_prefix_original_barcode_name_as_in_10x_directory
    doublet_cells <- new_cellnames[sel_index]
    print(head(doublet_cells))
    non_doublet_cells <- setdiff(new_cellnames,doublet_cells)
    doublet_cells_pred <- rep('yes',length(doublet_cells))
    names(doublet_cells_pred) <- doublet_cells
    non_doublet_cells_pred <- rep('no',length(non_doublet_cells))
    names(non_doublet_cells_pred) <- non_doublet_cells
    doublet_preds <- c(doublet_cells_pred,non_doublet_cells_pred)
    sel_cells <- intersect(Cells(seurat_object),names(doublet_preds))
    sel_doublet_preds <- doublet_preds[sel_cells]
    index <- match(Cells(seurat_object),names(sel_doublet_preds))
    mydoublet_pred_df[[sprintf("doublet_results_%s",i)]] <- sel_doublet_preds[index]
  }
  for(i in names(doublet_object[['doublet_scores']])){
    doubscores <- doublet_object[['doublet_scores']][[i]]
    names(doubscores) <- doublet_object$cellnames
    sel_cells <- intersect(Cells(seurat_object),names(doubscores))
    sel_doubscores <- doubscores[sel_cells]
    index <- match(Cells(seurat_object),names(sel_doubscores))
    mydoublet_score_df[[sprintf("doublet_scores_%s",i)]] <- sel_doubscores[index]
  }
  return(list(doublet_prediction=data.frame(mydoublet_pred_df),doublet_scores=data.frame(mydoublet_score_df)))
}


# df$count <- rowSums(df[c(1,3)] == "Yes")

input_df <- read.table(input_tsv_file,sep="\t",header=FALSE)
colnames(input_df) <- c('Sample','cellranger_10x_directory')

doublet_sample_list <- list()

for(i in 1:nrow(input_df)){
doublet_sample_list[[input_df$Sample[i]]] <- doublet_files[grep(input_df$Sample[i],doublet_files)]
}

seurat_object <- readRDS(Seurat_file)

seurat_object_list <- SplitObject(seurat_object, split.by = split_by)

doublet_prediction_df <- list()
doublet_scores_df <- list()


for(i in names(seurat_object_list)){
doublet_object <- readRDS(doublet_sample_list[[i]])
doublet_dfs <- add_doublet_predictions_to_seurat_singlesample_mod(seurat_object=seurat_object_list[[i]],doublet_object=doublet_object,seurat_cellname_prefix=i){
doublet_prediction_df[[i]] <- doublet_dfs$doublet_prediction
doublet_scores_df[[i]] <- doublet_dfs$doublet_scores
}


All_doublet_preds <- do.call(rbind,doublet_prediction_df)

ncol_doubs <- ncol(All_doublet_preds)
maj <- round((ncol_doubs/2)+1, digits = 0)

All_doublet_preds$prediction_aggrement <- rowSums(All_doublet_preds=='yes')
#doublet predictor aggrement for each cell

All_doublet_preds$majority_doublet_predictions <- 'no'
All_doublet_preds$majority_doublet_predictions[All_doublet_preds$prediction_aggrement >= maj]  <- 'yes'
#majority doublets : predictions for which majority of doublet predictors agree

All_doublet_preds$prediction_aggrement <- NULL

All_doublet_scores <- do.call(rbind,doublet_scores_df)

index <- match(Cells(seurat_object),rownames(All_doublet_preds))

All_doublet_preds <- All_doublet_preds[index,]

index <- match(Cells(seurat_object),rownames(All_doublet_scores))

All_doublet_scores[index,] <- All_doublet_scores[index,]

seurat_object <- AddMetaData(object = seurat_object,metadata =All_doublet_preds)

seurat_object <- AddMetaData(object = seurat_object,metadata =All_doublet_scores)

date = gsub(".rds",".doublet_calls.rds",basename(Seurat_file));

saveRDS(seurat_object, file = paste0('Cycling.SCT.PCA.UMAP.TSNE.CLUST.',".doublet_calls.",date,".rds"))

# test_doublet[['doublet_results']] <- FindDoublets(score.list=test_doublet$doublet_scores,rate=0.08)

# doublet_results_dfs <- add_doublet_predictions_to_seurat_singlesample_mod(seurat_object=OE1_plus_2,doublet_object=test_doublet,seurat_cellname_prefix='OE1_plus_2')
