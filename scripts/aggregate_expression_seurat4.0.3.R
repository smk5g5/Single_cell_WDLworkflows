.libPaths( c("/storage1/fs1/allegra.petti/Active/R_libs_scratch/RLibs_4.0.3",.libPaths()) )
library(Seurat)
library(dplyr)

get_avg_scaledexp <- function(seurat_object,clustering,sel_genes){
  DefaultAssay(seurat_object) <- 'RNA'
  top.genes <- unique(sel_genes)
  # specify clustering result to plot
  Idents(object=seurat_object) <- clustering; # make it the default identity - can be any Identity you choose
  clust.list=sort(levels(seurat_object)); # sort the clusters
  clust.id.field <- which(names(seurat_object@meta.data) == clustering);
  # extract the field that contains the cluster membership for each cell
  temp.means = list(); # initialize list of means for each cluster
  for (i in 1:length(clust.list)) {
    clust = clust.list[i];
    print (clust);
    cells.i <- as.factor(Cells(seurat_object)[which(seurat_object@meta.data[,clust.id.field] == clust)]); #
    seurat_object.i = FetchData(object=seurat_object, cells=cells.i, vars=top.genes, slot="data"); # use data slot from seurat_object. rows are cells, columns are genes
    temp.means[[i]] <- data.frame(colMeans(seurat_object.i))
  }
  clust.means = as.data.frame(dplyr::bind_cols(temp.means)); # rows are genes, columns are clusters
  names(clust.means) <- clust.list;
  rownames(clust.means) <- rownames(temp.means[[1]]);
  clust.means <- as.matrix(clust.means);
  rowmean = rowMeans(clust.means);
  rowsd = apply(clust.means, 1, sd);
  clust.means.norm = (clust.means-rowmean)/rowsd;
  return(clust.means.norm)
}


args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 1) {
  args <- c("--help")
}

Seurat_file <- as.character(args[1])
clustering <- as.character(args[2])
DEG_file <- as.character(args[3])


seurat_object <- readRDS(Seurat_file)

DEGs <- read.table(DEG_file,header=T)

DEGs_sig <- DEGs[DEGs$p_val_adj <= 0.05,]

top25_DEGs_sig <- DEGs_sig %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC))%>% 
  dplyr::slice(1:25)

clust.means.norm <- get_avg_scaledexp(seurat_object = seurat_object,clustering=clustering,sel_genes=unique(top25_DEGs_sig$gene))

saveRDS(clust.means.norm,'aggregate_expression_obj.rds')
