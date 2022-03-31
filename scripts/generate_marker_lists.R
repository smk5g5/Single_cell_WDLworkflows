df_to_marker_list <- function(marker_df,cluster_col,gene_col){
  sub_marker_list <- list()
  for(i in unique(marker_df[[cluster_col]])){
    subdf <- marker_df[marker_df[[cluster_col]]==i,]
    sub_marker_list[[i]]  <- subdf[[gene_col]]
  }
return(sub_marker_list)
  }

# get_topn_bycol <- function(mydf,group_col,arrange_col,slice_no){
# topn_mydf <- mydf %>%
#   group_by_at(group_col) %>%
#   arrange(desc(.dots=arrange_col)) %>%
#   slice(1:slice_no)
#   return(topn_mydf)
# }

Nowaski_DevelopingBrain <- read.table('/storage1/fs1/allegra.petti/Active/10xGenomics/key.gene.lists/brain_cell_markers/Nowaski_DevelopingBrain.txt',header=T)

Nowaski_DevelopingBrain_list <- df_to_marker_list(Nowaski_DevelopingBrain,cluster_col='List',gene_col='gene_name')

Li_etal <- read.table('/storage1/fs1/allegra.petti/Active/10xGenomics/key.gene.lists/brain_cell_markers/Li.txt',header=T,sep="\t")

Li_etal_list <- df_to_marker_list(Li_etal,cluster_col='List',gene_col='gene_name')

Ivy_markers <- read.table('/storage1/fs1/allegra.petti/Active/10xGenomics/key.gene.lists/brain_cell_markers/Ivy_markers.txt',header=T)

Ivy_markers_top25 <- Ivy_markers %>%
  group_by(List) %>%
  arrange(desc(minFoldChange)) %>%
  slice(1:25)

Ivy_markers_top25_list <- df_to_marker_list(Ivy_markers_top25,cluster_col='List',gene_col='gene_name')

Ivy_markers_top50 <- Ivy_markers %>%
  group_by(List) %>%
  arrange(desc(minFoldChange)) %>%
  slice(1:50)

Ivy_markers_top50_list <- df_to_marker_list(Ivy_markers_top50,cluster_col='List',gene_col='gene_name')

Wang_CancerCell <- read.table('/storage1/fs1/allegra.petti/Active/10xGenomics/key.gene.lists/brain_cell_markers/Wang_CancerCell.txt',header=T)

Wang_SingleCell_CancerDiscovery <- read.table('/storage1/fs1/allegra.petti/Active/10xGenomics/key.gene.lists/brain_cell_markers/Wang_SingleCell_CancerDiscovery2019.txt',header=T)

Wang_CancerCell_list <- df_to_marker_list(Wang_CancerCell,cluster_col='List',gene_col='gene_name')

Wang_SingleCell_CancerDiscovery_list <- df_to_marker_list(Wang_SingleCell_CancerDiscovery,cluster_col='List',gene_col='gene_name')

saveRDS(Ivy_markers_top25_list,'./essential_inputs/Ivy_markers_top25_list.rds')
saveRDS(Ivy_markers_top50_list,'./essential_inputs/Ivy_markers_top50_list.rds')
saveRDS(Li_etal_list,'./essential_inputs/Li_etal_list.rds')
saveRDS(Wang_CancerCell_list,'./essential_inputs/Wang_CancerCell_list.rds')
saveRDS(Wang_SingleCell_CancerDiscovery_list,'./essential_inputs/Wang_SingleCell_CancerDiscovery_list.rds')
saveRDS(Nowaski_DevelopingBrain_list,'./essential_inputs/Nowaski_DevelopingBrain_list.rds')


#/rdcw/fs1/allegra.petti/Active/10xGenomics/key.gene.lists/brain_cell_markers/Ivy_markers.txt

#/rdcw/fs1/allegra.petti/Active/10xGenomics/key.gene.lists/brain_cell_markers/Ivy_markers.txt

#

#/rdcw/fs1/allegra.petti/Active/10xGenomics/key.gene.lists/brain_cell_markers/Wang_SingleCell_CancerDiscovery2019.txt

############Clean Neftel data############

Neftel <- read.table('/storage1/fs1/allegra.petti/Active/10xGenomics/key.gene.lists/brain_cell_markers/neftel_gbm_signatures.txt',header=T,sep="\t")

clean_neftel_grep <- function(x){
return(x[grep('^$',x,invert=T)])
}

neftel_lists <- list()

for(i in colnames(Neftel)){
  neftel_lists[[i]] <- Neftel[[i]] 
}

neftel_lists_clnd <- lapply(neftel_lists,clean_neftel_grep)

saveRDS(neftel_lists_clnd,'./essential_inputs/Neftel_etal.rds')

###################################################

Panglioma <- read.table('/storage1/fs1/allegra.petti/Active/10xGenomics/key.gene.lists/brain_cell_markers/panglioma.txt',header=T,sep="\t")

panglioma_list <- df_to_marker_list(Panglioma,cluster_col='pan_glioma_state',gene_col='gene')

saveRDS(panglioma_list,'./essential_inputs/panglioma_list.rds')


mycelltype <- gsub(pattern = '\\/',replacement = '_or_',i)


gsub('.rds','',unlist(str_split(basename(outrds_name),pattern='_scsorter_preds_',simplify=T))[2])

