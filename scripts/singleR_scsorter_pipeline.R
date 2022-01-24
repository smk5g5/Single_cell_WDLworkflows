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
cluster.res =  as.numeric(args[3])
nPC= as.integer(args[4])
yaml_file =  as.character(args[5])

yaml_input <- args[1]

get_all_references <- function(yaml_obj){
  #test[['References']][[1]][['brain_im_atl']][['prefix']]
  all_refs <- list()
  for(i in 1:length(yaml_obj[['References']])){
    rn <- names(yaml_obj[['References']][[i]])
    pref <- yaml_obj[['References']][[i]][[rn]][['prefix']]
    fl  <- yaml_obj[['References']][[i]][[rn]][['file']]
    all_refs[[pref]] <- fl
  }
  ref_names <- names(all_refs)
  myref <- unlist(all_refs)
  return(myref)
}

get_sel_refs<- function(all_refs,yaml_obj,ref_type){
  sel_names <- yaml_obj[[ref_type]]
  sel_refs <- all_refs[sel_names]
  return(sel_refs)
}

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


date = gsub("2021-","21",Sys.Date(),perl=TRUE);
date = gsub("-","",date);

seurat_obj <- readRDS(Seurat_rds)



cells_by_cluster_percent <- function(seurat_obj,cluster_meta,sel_cells,cell_type,kmeans_1d=TRUE){
  Idents(seurat_obj) <- cluster_meta
  cluster_ids <- names(table(Idents(seurat_obj)))
  Cells_by_cluster <- data.frame(Cluster=character(0),fraction=numeric(0))
  counter=0
  for(i in cluster_ids){
    seurat_obj_sub <- subset(seurat_obj,idents=i)
    common_cells <- intersect(sel_cells,Cells(seurat_obj_sub))
    cell_frac <- length(common_cells)/length(Cells(seurat_obj_sub))
    print(cell_frac)
    print(i)
    counter= counter+1
    Cells_by_cluster[counter, ] <- c(i,cell_frac)
  }
  Cells_by_cluster$fraction <- as.numeric(Cells_by_cluster$fraction)
  if(kmeans_1d=='TRUE'){
    #if kmeans_1d is True do univariate clustering of the fraction column with K=2
    # result <- Ckmeans.1d.dp(x = myexp$RNA,k = 2)
    result <- Ckmeans.1d.dp(x = Cells_by_cluster$fraction,k = 2)
    Cells_by_cluster$kmeans_cluster <- result$cluster
    if(max(Cells_by_cluster$fraction[Cells_by_cluster$kmeans_cluster==1]) > max(Cells_by_cluster$fraction[Cells_by_cluster$kmeans_cluster==2])){
      Cells_by_cluster$Celltype[Cells_by_cluster$kmeans_cluster==1] <- cell_type
      Cells_by_cluster$Celltype[Cells_by_cluster$kmeans_cluster==2] <- 'No.prediction'
    }
    else{
      Cells_by_cluster$Celltype[Cells_by_cluster$kmeans_cluster==2] <- cell_type
      Cells_by_cluster$Celltype[Cells_by_cluster$kmeans_cluster==1] <- 'No.prediction'
    }
    return(Cells_by_cluster)
  }
  else{
  max_val <- max(mean(Cells_by_cluster$fraction),median(Cells_by_cluster$fraction))
  Cells_by_cluster$Celltype[Cells_by_cluster$fraction >= max_val] <- cell_type
  Cells_by_cluster$Celltype[Cells_by_cluster$fraction <= max_val] <- 'No.prediction'
  return(Cells_by_cluster)
  }
}

jpeg(filename="B189_cd45pos_clusters.jpg", units="in", res=300, height=11, width=14);
B189_cd45pos <- DimPlot(B189, cells.highlight= Cells(subset(B189,idents=immune_clus)), cols.highlight = 'red',cols="grey88",sizes.highlight=0.1,order=T,label=T,label.size=5,raster=F)+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank())+ggtitle('CD45+ve clusters') 
print(B189_cd45pos)
dev.off()

jpeg(filename="B189_cd45neg_clusters.jpg", units="in", res=300, height=11, width=14);
B189_cd45neg <- DimPlot(B189, cells.highlight= Cells(subset(B189,idents=nonimmune_clus)), cols.highlight = 'red',cols="grey88",sizes.highlight=0.1,order=T,label=T,label.size=5,raster=F)+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank())+ggtitle('CD45-ve clusters')
print(B189_cd45neg)
dev.off()



