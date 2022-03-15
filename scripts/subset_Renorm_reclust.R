.libPaths( c("/storage1/fs1/allegra.petti/Active/R_libs_scratch/RLibs_4.0.3",.libPaths()) )
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(Matrix)
library(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 5) {
  args <- c("--help")
}

seurat_loc <- as.character(args[1])
sub_col <- as.character(args[2])
tech <- as.character(args[3])
inverse <- as.character(args[3])
output.stats <- as.character(args[4])
output_meta <- as.character(args[5])

seurat_obj <- readRDS(seurat_loc)

print(names(seurat_obj@meta.data))

subset_renormalize_recluster<- function(seurat_obj,sub_col,tech,inverse,output_meta) {
  DefaultAssay(seurat_obj) <- "RNA"
  Idents(seurat_obj) <- sub_col
  if(inverse==TRUE){
    scrna_GEX <- subset(seurat_obj,idents=tech,invert = TRUE)
  }else{
    scrna_GEX <- subset(seurat_obj,idents=tech)
  }
  scrna_GEX <- FindVariableFeatures(scrna_GEX)  
  genome <- "GRCh38";
  date = gsub("2022-","22",Sys.Date(),perl=TRUE);
  date = gsub("-","",date);

  scrna_GEX <- ScaleData(scrna_GEX, verbose = FALSE)
  scrna_GEX <- RunPCA(scrna_GEX, npcs = 30, verbose = FALSE)
  scrna_GEX <- RunUMAP(scrna_GEX, reduction = "pca", dims = 1:30)
  scrna_GEX <- FindNeighbors(scrna_GEX, reduction = "pca", dims = 1:30)
  scrna_GEX <- FindClusters(scrna_GEX, resolution = c(0.5, 0.7, 0.9))
  
  scrna_GEX <- StashIdent(object = scrna_GEX, save.name = sprintf("Clusters_scrnaonly_%.1f_%dPC",0.5, 30))
  scrna_GEX <- StashIdent(object = scrna_GEX, save.name = sprintf("Clusters_scrnaonly_%.1f_%dPC",0.7, 30))
  scrna_GEX <- StashIdent(object = scrna_GEX, save.name = sprintf("Clusters_scrnaonly_%.1f_%dPC",0.9, 30))
  Idents(scrna_GEX) <- sprintf("Clusters_scrnaonly_%.1f_%dPC",0.5, 30)
  DEGs <- FindAllMarkers(object=scrna_GEX); # output is a matrix!
  write.table(DEGs, file=sprintf("%s/DEGs.Wilcox.PCA.%d.cluster.%.1f.%s.xls", output.stats, 30, 0.5, date), quote=FALSE, sep="\t", row.names=FALSE) # must save cluster-specific marker genes
	
  n.graph = length(unique(Idents(scrna_GEX)))
  print(n.graph)
  rainbow.colors = rainbow(n.graph, s=0.6, v=0.9);
  names(rainbow.colors) <- sort(unique(Idents(scrna_GEX)))
  print(rainbow.colors)
  print(length(rainbow.colors))

  jpeg(sprintf("%s/UMAP.clusters.%d.%.1f.%s.jpg", output.stats, 30, 0.5, date), width = 10, height = 8, units="in", res=300);
  p2 <- DimPlot(object = scrna_GEX, reduction = "umap", group.by = sprintf("Clusters_scrnaonly_%.1f_%dPC",0.5, 30), cols = rainbow.colors, pt.size=0.1) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  print(p2);
  dev.off();
  
  jpeg(sprintf("%s/UMAP.clusters.labeled.%d.%.1f.%s.jpg", output.stats,30, 0.5, date), width = 10, height = 8, units="in", res=300);
  p2 <- DimPlot(object = scrna_GEX, reduction = "umap", group.by = sprintf("Clusters_scrnaonly_%.1f_%dPC",0.5, 30), cols = rainbow.colors, pt.size=0.1, label=TRUE,label.size = 5) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  print(p2);
  dev.off();
  
  print ("color UMAP by Principal Components");
  jpeg(sprintf("%s/UMAP.%d.%.1f.colorby.PCs.%s.jpg", output.stats, 30, 0.5, date), width = 12, height = 6, units="in", res=100);
  redblue=c("blue","gray","red");
  fp1 <- FeaturePlot(object = scrna_GEX, features = 'PC_1', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  fp2 <- FeaturePlot(object = scrna_GEX, features = 'PC_2', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  fp3 <- FeaturePlot(object = scrna_GEX, features = 'PC_3', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  fp4 <- FeaturePlot(object = scrna_GEX, features = 'PC_4', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  fp5 <- FeaturePlot(object = scrna_GEX, features = 'PC_5', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  fp6 <- FeaturePlot(object = scrna_GEX, features = 'PC_6', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  fp7 <- FeaturePlot(object = scrna_GEX, features = 'PC_7', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  fp8 <- FeaturePlot(object = scrna_GEX, features = 'PC_8', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  fp9 <- FeaturePlot(object = scrna_GEX, features = 'PC_9', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  fp10 <- FeaturePlot(object = scrna_GEX, features = 'PC_10', cols=redblue, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  print(plot_grid(fp1, fp2, fp3, fp4, fp5, fp6, fp7, fp8, fp9, fp10));
  # print(plot_grid(fp1, fp2, fp3, fp4, fp5, fp6, fp7, fp8, fp9));
  dev.off();
  output_name <- paste0(output.stats,'/',output_meta,'.RDS')
  saveRDS(scrna_GEX,file=output_name)
}

subset_recluster_renormalize(seurat_obj=seurat_obj,sub_col=sub_col,tech=tech,output.stats=output.stats,output_meta=output_meta)
