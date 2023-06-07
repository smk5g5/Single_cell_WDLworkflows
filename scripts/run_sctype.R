.libPaths( c("/storage1/fs1/allegra.petti/Active/R_libs_scratch/RLibs_4.0.3",.libPaths()) )
library(Seurat)
library(BiocParallel)
library(yaml)
library(Ckmeans.1d.dp)
library(stringr)
library(scSorter)
require("biomaRt")
library(RColorBrewer)
library(dplyr)
library(ggplot2)

date = gsub("2022-","22",Sys.Date(),perl=TRUE);
date = gsub("-","",date);

args = commandArgs(trailingOnly=TRUE);
print (length(args))

if(length(args) < 4)
{
  stop("\n Number of inputs provided not enough \n")
}

seurat_rds <- args[1]
cluster_var <- args[2]
sctype_meta_name <- args[3]
output_name <- args[4]

seurat_obj <- readRDS(seurat_rds)

lapply(c("dplyr","Seurat","HGNChelper"), library, character.only = T)

run_sctype <- function(seurat_obj,cluster_var,sctype_meta_name){
  Idents(seurat_obj) <- cluster_var
  #brain and immune together
  db_ = "https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_full.xlsx";
  gs_list_combined <- list()
  gs_list_immune = gene_sets_prepare(db_,"Immune system")
  gs_list_brain = gene_sets_prepare(db_, 'Brain')
  print('getting to making gs_list_combined')
  gs_list_combined[['gs_positive']] = c(gs_list_immune$gs_positive,gs_list_brain$gs_positive)
  gs_list_combined[['gs_negative']] = c(gs_list_immune$gs_negative,gs_list_brain$gs_negative)
  print('getting to using gs_list_combined')
  es.max_immune = sctype_score(scRNAseqData = seurat_obj[["RNA"]]@data, scaled = F, 
                               gs = gs_list_combined$gs_positive, gs2 = gs_list_combined$gs_negative)
  cL = do.call("rbind", lapply(as.character(unique(seurat_obj@meta.data[[cluster_var]])), function(cl){
    es.max.cl = sort(rowSums(es.max_immune[ ,rownames(seurat_obj@meta.data[seurat_obj@meta.data[cluster_var]==cl, ])]), decreasing = !0)
    head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_obj@meta.data[cluster_var]==cl)), 10)
  }))
    #   cL = do.call("rbind", lapply(unique(seurat_obj@meta.data[[cluster_var]]), function(cl){
  #   es.max.cl = sort(rowSums(es.max_immune[ ,rownames(seurat_obj@meta.data[seurat_obj@meta.data[cluster_var]==cl, ])]), decreasing = !0)
  #   head(data.frame(cluster = cl, type = names(es.max.cl), scores = es.max.cl, ncells = sum(seurat_obj@meta.data[cluster_var]==cl)), 10)
  # }))
  # 
  sctype_scores_immune  = cL %>% group_by(cluster) %>% top_n(n = 1, wt = scores)  
  
  sctype_scores_immune$type[as.numeric(as.character(sctype_scores_immune$scores)) < sctype_scores_immune$ncells/4] = "Unknown"
  
  seurat_obj@meta.data[sctype_meta_name] = ""
  for(j in unique(sctype_scores_immune$cluster)){
    cl_type = sctype_scores_immune[sctype_scores_immune$cluster==j,]; 
    seurat_obj@meta.data[sctype_meta_name][seurat_obj@meta.data[cluster_var] == j] = as.character(cl_type$type[1])
  }
  
  return(list('seurat'=seurat_obj,'sctype_res'=cL,'sctype_scores'=sctype_scores_immune))
}

sctype_res_list <- run_sctype(seurat_obj=seurat_obj,cluster_var=cluster_var,sctype_meta_name=sctype_meta_name)

saveRDS(sctype_res_list[['seurat']],file=paste0(output_name,'_sctype_results.',date,'.rds'))
