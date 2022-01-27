.libPaths( c("/storage1/fs1/allegra.petti/Active/R_libs_scratch/RLibs_4.0.3",.libPaths()) )
library(clusterProfiler)
library(org.Hs.eg.db)


make_logfc_sorted_genelist_bycluster <- function(DEG_results){
  sorted_fc_list <- list()
  for(i in unique(DEG_results$cluster)){
    clus_df <- subset(DEG_results,cluster==i)
    gene_vec <- clus_df$avg_log2FC
    names(gene_vec) <- clus_df$gene
    gene_vec<-na.omit(gene_vec)
    gene_vec = sort(gene_vec, decreasing = TRUE)
    sorted_fc_list[[i]] <- gene_vec
  }
  return(sorted_fc_list)
}
#DEGs is the FindAllmarkers result object from Seurat
fc_list <- make_logfc_sorted_genelist_bycluster(DEGs)

GSEA_run_func <- function(x,myorg){
  gse <- gseGO(geneList=x, 
               ont ="BP", 
               keyType = "SYMBOL", 
               minGSSize = 5, 
               maxGSSize = 500, 
               pvalueCutoff = 0.05,
               pAdjustMethod = "BH",
               verbose = TRUE, 
               OrgDb = myorg)
  return(gse)
}

#Run GSEA analysis on each cluster

gse_list_byclus <- lapply(fc_list,GSEA_run_func,myorg="org.Hs.eg.db")

names(gse_list_byclus) <- paste0('cluster_',names(gse_list_byclus))

#Merge gsea results into one object
gsea_merged_results <- simplify(merge_result(gse_list_byclus_simp))

#write GSEA results to outfile
write.table(gsea_merged_results@compareClusterResult, file = 'Clusterwide_GSEA_GOBP.txt',
            row.names = F,sep = "\t",quote = F) 

#enrichment of msigdb terms
#function takes three arguments x is the gene list cat is the category of msigdb and 
#sp is the species e.g. "Homo sapiens"

msigdb_enrichment <- function(x,cat,sp){
  if(cat=="C5"){
    m_t2g <- msigdbr(species = sp, category = "C5",subcategory = "GO:BP") %>% 
      dplyr::select(gs_name, gene_symbol)
  }
  else if(cat=="C3"){
    m_t2g <- rbind(msigdbr(species = sp, category = "C3",subcategory = 'TFT:TFT_Legacy'),msigdbr(species = sp, category = "C3",subcategory = 'TFT:GTRD'))  %>% 
      dplyr::select(gs_name, gene_symbol)
  }
  else if(cat=="H"){
    m_t2g <- msigdbr(species = sp, category = "H") %>% 
      dplyr::select(gs_name, gene_symbol)
  }
  else if(cat=="C6"){
    m_t2g <- msigdbr(species = sp, category = "C6") %>% 
      dplyr::select(gs_name, gene_symbol)
  }
  else if(cat=="C7"){
    m_t2g <- msigdbr(species = sp, category = "C7") %>% 
      dplyr::select(gs_name, gene_symbol)
  }
  else{
    m_t2g <- msigdbr(species = sp, category = "C8") %>% 
      dplyr::select(gs_name, gene_symbol)
  }
  em2 <- GSEA(x, TERM2GENE = m_t2g,minGSSize = 10,maxGSSize = 500, pvalueCutoff =1)
  return(em2)
}
