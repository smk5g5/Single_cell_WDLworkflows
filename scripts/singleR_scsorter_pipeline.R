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

################################
##This script reads in the parameters from yaml file and divides the data first into ptprc+ve/-ve
##clusters then uses immune references to make predictions using SingleR and non immune references
##to make predictions using SingleR
##################################

#Naming convention of cluster metadata column
#ClusterNames_0.3_15PC
#sprintf("ClusterNames_%.1f_%dPC", cluster.res, nPC)

######################################
#further subdivide malignant cells into subtypes it takes a
##list of markers by subtypes and make predictions on malignant 
#cell subtypes using scsorter more details about the parameters 
#are in the yaml file
################################

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 5) {
  args <- c("--help")
}

Seurat_rds <- as.character(args[1])
sample_name <- as.character(args[2])
cluster.res =  as.numeric(args[3])
nPC= as.integer(args[4])
yaml_input =  as.character(args[5])

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

cells_by_cluster_percent <- function(seurat_obj,cluster_meta,sel_cells,cell_type,imm_params){
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
  if(imm_params['kmeans_1d']=='TRUE'){
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

make_singleR_refs <- function(ref_vec){
  singleR_ref_list <- list()
 for(i in names(ref_vec)){
   singleR_ref_list[[i]] <- readRDS(ref_vec[i])
   singleR_ref_list[[i]]$label <- paste(i,singleR_ref_list[[i]]$label,sep='.')
 }
  return(singleR_ref_list)
}

make_singleR_labels <- function(singleR_ref_list){
  label_list <- list()
  for(i in 1:length(names(singleR_ref_list))){
    ref_name <- names(singleR_ref_list)[i]
    label_list[[i]] <- singleR_ref_list[[ref_name]]$label
  }
  return(label_list)
}


date = gsub("2021-","21",Sys.Date(),perl=TRUE);
date = gsub("-","",date);

seurat_obj <- readRDS(Seurat_rds)

Cluster_meta_name <- sprintf("ClusterNames_%.1f_%dPC", cluster.res, nPC)

#read Yaml object
yaml_obj <- read_yaml(yaml_input)
#parse references from yaml object
all_refs <- get_all_references(yaml_obj= yaml_obj)

#parse Immune_reference_names from yaml object
immune_refs <- get_sel_refs(all_refs = all_refs,yaml_obj=yaml_obj,ref_type='Immune_reference_names')

#parse NonImmune_reference_names from yaml object
non_immune_refs <- get_sel_refs(all_refs = all_refs,yaml_obj=yaml_obj,ref_type='NonImmune_reference_names')

#parse Immune_Cell_Selection parameters from yaml object
imm_params <- unlist(yaml_obj[['Immune_Cell_Selection']])

ptprc_counts <- as.matrix(seurat_obj@assays$RNA@counts)['PTPRC',]
ptprc_cells <- names(ptprc_counts[ptprc_counts >= 1])

#output fraction of cluster percentages containg ptprc+ve cells
#cluster_df has 4 columns (Cluster,fraction,kmeans_cluster,Celltype) if kmeans_1d is True otherwise 3 columns (Cluster,fraction,Celltype)

cluster_df <- cells_by_cluster_percent(seurat_obj=seurat_obj,cluster_meta=Cluster_meta_name,sel_cells=ptprc_cells,cell_type=yaml_obj[['Subset_celltype']],imm_params=imm_params)

#select clusters from seurat object which are immune cell clusters

sel_clus <- cluster_df$Cluster[cluster_df$Celltype==yaml_obj[['Subset_celltype']]]

#make the Cluster_meta_name as the new ident for seurat object
Idents(seurat_obj) <- Cluster_meta_name

jpeg(filename=sprintf("%s_cd45pos_clusters.jpg", sample_name), units="cm", res=300, height=11, width=14);
cd45pos_clus <- DimPlot(seurat_obj, cells.highlight= Cells(subset(seurat_obj,idents=sel_clus)), cols.highlight = 'red',cols="grey88",sizes.highlight=0.1,order=T,label=T,label.size=5,raster=F)+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank())+ggtitle('CD45+ve clusters') 
print(cd45pos_clus)
dev.off()

jpeg(filename=sprintf("%s_cd45neg_clusters.jpg", sample_name), units="cm", res=300, height=11, width=14);
cd45neg_clus <- DimPlot(seurat_obj, cells.highlight= Cells(subset(seurat_obj,idents=sel_clus,invert=T)), cols.highlight = 'red',cols="grey88",sizes.highlight=0.1,order=T,label=T,label.size=5,raster=F)+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank())+ggtitle('CD45-ve clusters')
print(cd45neg_clus)
dev.off()

#select immune clusters
seurat_obj_immune <- subset(seurat_obj,idents=sel_clus)

#read the immune singleR references into a list
immune_ref_list <- make_singleR_refs(immune_refs)

#make SingleR predictions on immune cell clusters
singleR_immune_res <- SingleR(test = as.SingleCellExperiment(seurat_obj_immune),
                             ref = immune_ref_list,
                             labels = make_singleR_labels(immune_ref_list),
                             BPPARAM=MulticoreParam())

saveRDS(object=singleR_immune_res,file=sprintf("%s_singleR_immune_res.rds", sample_name))

#saveRDS(abc_momac_verse_diff,'mo_mac_verse_mulderetal_comp.rds')

#select non immune clusters
seurat_obj_nonimmune <- subset(seurat_obj,idents=sel_clus,invert=T)

#read the nonimmune singleR references into a list
nonimmune_ref_list <- make_singleR_refs(non_immune_refs)

if(yaml_obj[['Merge_Neftel_refs']]=='TRUE'){
  for(i in names(non_immune_refs[grep('nef',names(non_immune_refs))])){
    nonimmune_ref_list[[i]]= nonimmune_ref_list[[i]][,grep('Malignant',colnames(nonimmune_ref_list[[i]]))]
    nonimmune_ref_list[[i]]$label <- str_replace(nonimmune_ref_list[[i]]$label,i,yaml_obj[['Malignant_prefix']])
  }
}

#Run SingleR on non immune subset of cells
singleR_nonimmune_res <- SingleR(test = as.SingleCellExperiment(seurat_obj_nonimmune),ref = nonimmune_ref_list,
                                labels = make_singleR_labels(nonimmune_ref_list),BPPARAM=MulticoreParam())

saveRDS(object=singleR_nonimmune_res,file=sprintf("%s_singleR_nonimmune_res.rds", sample_name))