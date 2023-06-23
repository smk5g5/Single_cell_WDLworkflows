library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(Matrix)
library("sctransform");
library("dplyr");
library("RColorBrewer");
library("ggplot2");
library("cowplot");
library(tidyverse)

library(ggplot2)
library(gridExtra)
library(ggthemes)


args <- commandArgs(trailingOnly = TRUE)

project_name <- as.character(args[1])
organism <- as.character(args[2])
seurat_files <- as.character(args[3:length(args)])

if(organism=='human'){
genome <- "GRCh38";
} else{
  genome <- "GRCm38";
}
date = gsub("2023-","23",Sys.Date(),perl=TRUE);
date = gsub("-","",date);
# 
# recale_plot_agg <- function(myplot_obj,plot_name,mywid,myht,mydpi) {
#   ggsave(plot_name, plot = myplot_obj, device = ragg::agg_png, width = mywid, height = myht, units = "in",dpi = mydpi) # units = "in"
#   ggs <- image_read(plot_name)
#   print(ggs)
# }


options(future.globals.maxSize = 5000 * 1024^2)
plan("multiprocess", workers = 4)


get_significant_pcs <- function(scrna_GEX) {
  control='Cycling'
  # if(length(unique(scrna_GEX$orig.ident))==1){
  #  nPC=20 
  # }else{
    nPC=100
  # }
  # print('Does the error happen inside get_significant_pcs?')
  scrna_GEX <- RunPCA(scrna_GEX, npcs = nPC, verbose = FALSE)
  scrna_GEX <- JackStraw(object = scrna_GEX, num.replicate = 100, dims=nPC)
  scrna_GEX <- ScoreJackStraw(object = scrna_GEX, dims = 1:nPC)
  jpeg(sprintf("PCA.jackstraw.%s.%s.%s.jpg", control,project_name,date), width = 10, height = 6, units="in", res=300);
  js <- JackStrawPlot(object = scrna_GEX, dims = 1:nPC)
  print(js);
  dev.off();
  pc.pval <- scrna_GEX@reductions$pca@jackstraw@overall.p.values
  print(pc.pval);
  nPC=length( pc.pval[,'Score'][pc.pval[,'Score'] <= 0.05]) 
  # print('No the error does not happen inside get_significant_pcs!')
  #redefine nPCs based on number of significant prinicipal components in jackstraw plot
  return(nPC)
}

renormalize_recluster <- function(scrna_GEX,date,organism) {
  # DefaultAssay(seurat_obj) <- "RNA"
  # Idents(seurat_obj) <- sub_col
  # if(inverse=='TRUE'){
  #   scrna_GEX <- subset(seurat_obj,idents=ident_names,invert = TRUE)
  # }else{
  #   scrna_GEX <- subset(seurat_obj,idents=ident_names)
  # }

  # ##################################################################################
  # # This bit is only for testing purposes for the rscript within wdl would be disabled
  # # in the main workflow
  # set.seed(100)
  # random_sample_of_cells = sample(Cells(scrna_GEX),length(Cells(scrna_GEX)) * 0.1)
  # #select 10% of all cells randomly for testing the script.
  # scrna_GEX <- subset(scrna_GEX,cells=random_sample_of_cells)
  # ##################################################################################

  if(organism=='human'){
  cell.cycle.tirosh <- read.table("/n/scratch3/users/s/sak4832/meningioma_project/For_cellbender/CellCycleTirosh.txt", sep='\t', header=TRUE);
  s.genes = cell.cycle.tirosh$`Gene.Symbol`[which(cell.cycle.tirosh$List == "G1/S")];
  g2m.genes = cell.cycle.tirosh$`Gene.Symbol`[which(cell.cycle.tirosh$List == "G2/M")];
  } else{
  cell.cycle.tirosh <- read.table("/n/scratch3/users/s/sak4832/meningioma_project/For_cellbender/CellCycleTirosh_mouse.txt", sep='\t', header=FALSE);
  s.genes = cell.cycle.tirosh$V2[which(cell.cycle.tirosh$V1 == "G1/S")];
  g2m.genes = cell.cycle.tirosh$V2[which(cell.cycle.tirosh$V1 == "G2/M")];
  }

  scrna_GEX <- CellCycleScoring(object=scrna_GEX, s.features=s.genes, g2m.features=g2m.genes, set.ident=FALSE)

  scrna_GEX <- NormalizeData(object = scrna_GEX, normalization.method = "LogNormalize", scale.factor = 1e6); # 1e6 is new as of 1/8/20

  scrna_GEX <- FindVariableFeatures(object = scrna_GEX, selection.method = 'vst', mean.cutoff = c(0.1,8), dispersion.cutoff = c(1, Inf))

  control='Cycling'

  if (control == "Cycling") { # This removes all signal associated with the cell cycle
   scrna_GEX <- ScaleData(object = scrna_GEX, features = rownames(x = scrna_GEX), vars.to.regress = c("S.Score","G2M.Score"), display.progress=FALSE);
  } else if (control == "CyclingRB") {
   scrna_GEX <- ScaleData(object = scrna_GEX, features = rownames(x = scrna_GEX), vars.to.regress = c("S.Score","G2M.Score","percent.ribo"), display.progress=FALSE);
  } else if (control == "CyclingDiff") {
   scrna_GEX <- ScaleData(object = scrna_GEX, features = rownames(x = scrna_GEX), vars.to.regress = c("CC.Difference"), display.progress=FALSE);
  } else {
   scrna_GEX <- ScaleData(object = scrna_GEX, features = rownames(x = scrna_GEX), display.progress=FALSE);
  }

  nPC <- get_significant_pcs(scrna_GEX)
  
  scrna_GEX <- RunPCA(scrna_GEX, npcs = nPC, verbose = FALSE)

  scrna_GEX <- RunUMAP(object = scrna_GEX, reduction = "pca", dims = 1:nPC)
  scrna_GEX <- RunTSNE(object = scrna_GEX, reduction = "pca", dims = 1:nPC)

  scrna_GEX <- JackStraw(object = scrna_GEX, num.replicate = 100, dims=nPC)
  scrna_GEX <- ScoreJackStraw(object = scrna_GEX, dims = 1:nPC)
  jpeg(sprintf("PCA.jackstraw.%s.%s.%s.jpg",control,project_name,date), width = 10, height = 6, units="in", res=300);
  js <- JackStrawPlot(object = scrna_GEX, dims = 1:nPC)
  print(js);
  dev.off();

  jpeg(sprintf("UMAP.%s.%s.%s.jpg",control,project_name, date), width = 10, height = 8, units="in", res=300);
  p2 <- DimPlot(object = scrna_GEX, reduction = "umap", group.by = "Sample", pt.size=0.1)
  print(p2);
  dev.off();

  print ("VizDimLoadings Running...");
  jpeg(sprintf("VizDimLoadings.%s.%s.%s.jpg",control,project_name,date), width = 8, height = 30, units="in", res=300);
  vdl <- VizDimLoadings(object = scrna_GEX, dims = 1:3)
  print(vdl);
  dev.off();

  print ("ProjectDim Running...");
  scrna_GEX <- ProjectDim(object = scrna_GEX)

  # saveRDS(scrna_GEX, file = sprintf("%s.SCT.PCA.UMAP.TSNE.%s.rds",control, date))

  print ("DimHeatmap Running...");
  jpeg(sprintf("PCA.heatmap.top.%s.%s.%s.jpg",control,project_name,date), width = 8.5, height = 11, units="in", res=300);
  hm <- DimHeatmap(object = scrna_GEX, dims = 1, cells = 500, balanced = TRUE);
  print(hm);
  dev.off();

  jpeg(sprintf("PCA.heatmap.multi.%s.%s.%s.jpg",control,project_name,date), width = 8.5, height = 24, units="in", res=300);
  hm.multi <- DimHeatmap(object = scrna_GEX, dims = 1:10, cells = 500, balanced = TRUE);
  print(hm.multi);
  dev.off();

  scrna_GEX <- FindNeighbors(scrna_GEX, reduction = "pca", dims = 1:nPC)
  scrna_GEX <- FindClusters(scrna_GEX, resolution = 0.5)
  scrna_GEX[[sprintf("ClusterNames_%.1f_%dPC",0.5, nPC)]] <- Idents(scrna_GEX)
  scrna_GEX <- FindClusters(scrna_GEX, resolution = 0.7)
  scrna_GEX[[sprintf("ClusterNames_%.1f_%dPC",0.7, nPC)]] <- Idents(scrna_GEX)
  scrna_GEX <- FindClusters(scrna_GEX, resolution = 0.9)
  scrna_GEX[[sprintf("ClusterNames_%.1f_%dPC",0.9, nPC)]] <- Idents(scrna_GEX)

  scrna_GEX <- FindClusters(scrna_GEX, resolution = 1.2)
  scrna_GEX[[sprintf("ClusterNames_%.1f_%dPC",1.2, nPC)]] <- Idents(scrna_GEX)

  #for now 0.7 is used as a default cluster resolution but something 
  #better needs to be replace this. Perhaps will use multiK 
  #https://www.ncbi.nlm.nih.gov/pmc/articles/PMC8375188/
  #https://github.com/siyao-liu/MultiK
  #in the next version 

  # Idents(scrna_GEX) <- sprintf("ClusterNames_%.1f_%dPC",1.2, nPC)
  # DEGs <- FindAllMarkers(object=scrna_GEX); # output is a matrix!
  # write.table(DEGs, file=sprintf("DEGs.Wilcox.PCA.%d.cluster.%.1f.%s.xls", nPC, 1.2, date), quote=FALSE, sep="\t", row.names=FALSE) # must save cluster-specific marker genes
  # 
  n.graph = length(unique(Idents(scrna_GEX)))
  print(n.graph)
  rainbow.colors = rainbow(n.graph, s=0.6, v=0.9);
  names(rainbow.colors) <- sort(unique(Idents(scrna_GEX)))
  print(rainbow.colors)
  print(length(rainbow.colors))

  jpeg(sprintf("UMAP.clusters.%d.%.1f.%s.%s.jpg", nPC, 1.2,project_name,date), width = 10, height = 8, units="in", res=300);
  p2 <- DimPlot(object = scrna_GEX, reduction = "umap", group.by = sprintf("ClusterNames_%.1f_%dPC",0.7, nPC), cols = rainbow.colors, pt.size=0.1) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  print(p2);
  dev.off();
  
  jpeg(sprintf("UMAP.clusters.labeled.%d.%.1f.%s.%s.jpg",nPC, 0.7,project_name,date), width = 10, height = 8, units="in", res=300);
  p2 <- DimPlot(object = scrna_GEX, reduction = "umap", group.by = sprintf("ClusterNames_%.1f_%dPC",0.7, nPC), cols = rainbow.colors, pt.size=0.1, label=TRUE,label.size = 5) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  print(p2);
  dev.off();

  feature.pal = rev(colorRampPalette(brewer.pal(11,"Spectral"))(50));
  
  mito.genes <- grep(pattern = "^MT-", x = rownames(x = scrna_GEX), value = TRUE,ignore.case = TRUE);
  percent.mito <- Matrix::colSums(x = GetAssayData(object = scrna_GEX, slot = 'counts')[mito.genes, ]) / Matrix::colSums(x = GetAssayData(object = scrna_GEX, slot = 'counts'));
  scrna_GEX[['percent.mito']] <- percent.mito;

  # ribosomal genes
  ribo.genes <- grep(pattern = "^RP[SL][[:digit:]]", x = rownames(x = scrna_GEX), value = TRUE,ignore.case = TRUE);
  percent.ribo <- Matrix::colSums(x = GetAssayData(object = scrna_GEX, slot = 'counts')[ribo.genes, ]) / Matrix::colSums(x = GetAssayData(object = scrna_GEX, slot = 'counts'));
  scrna_GEX[['percent.ribo']] <- percent.ribo;

  scrna_GEX[['percent.ribo']] <- percent.ribo;


  print("Making additional UMAP plots");
  # color UMAP plots by parameters of interest:
  print ("color by UMI");
  jpeg(sprintf("umap.%d.%.1f.colorby.UMI.%s.%s.%s.jpg", nPC, 0.7, control,project_name,date), width = 10, height = 8,  units="in", res=300);
  fp2 <- FeaturePlot(object = scrna_GEX, features = c("nCount_RNA"), cols = feature.pal, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  print(fp2);
  dev.off();

  print ("color by % mito");
  jpeg(sprintf("umap.%d.%.1f.colorby.MC.%s.%s.%s.jpg", nPC, 0.7, control,project_name,date), width = 10, height = 8, units="in", res=300);
  fp2 <- FeaturePlot(object = scrna_GEX, features = c("percent.mito"), cols = feature.pal, pt.size=0.1, reduction = "umap")+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  print(fp2);
  dev.off();

  print ("color by % RB");
  jpeg(sprintf("umap.%d.%.1f.colorby.RB.%s.%s.%s.jpg", nPC, 0.7, control,project_name,date), width = 10, height = 8, units="in", res=300);
  fp2 <- FeaturePlot(object = scrna_GEX, features = c("percent.ribo"), cols = feature.pal, pt.size=0.1, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  print(fp2);
  dev.off();

  print ("color by nGene");
  jpeg(sprintf("umap.%d.%.1f.colorby.nGene.%s.%s.%s.jpg", nPC, 0.7, control,project_name,date), width = 10, height = 8, units="in", res=300);
  fp2 <- FeaturePlot(object = scrna_GEX, features = c("nFeature_RNA"), cols = feature.pal, pt.size=0.1, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  print(fp2);
  dev.off();

  print ("color by Phase");
  jpeg(sprintf("umap.%d.%.1f.colorby.Phase.%s.%s.%s.jpg", nPC, 0.7, control,project_name,date,), width = 10, height = 8, units="in", res=300);
  phase.colors = ptol_pal()(3)
  umapplot <- DimPlot(object = scrna_GEX, cols=phase.colors, group.by="Phase", pt.size=0.1, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank())
  print(umapplot);
  dev.off();

  
  print ("color UMAP by Principal Components");
  jpeg(sprintf("UMAP.%d.%.1f.colorby.PCs.%s.%s.jpg", nPC, 0.7,project_name,date), width = 12, height = 6, units="in", res=100);
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
  # output_name <- paste0(output.stats,'/',output_meta,'.RDS')
  return(scrna_GEX)
}

seurat_obj_list <- list()
sample_names <- c()

# mysrt_files = read.table(seurat_files,header = F)

for(i in seurat_files){
  obj_name <- i
  test_obj <- readRDS(obj_name)
  sample_names[i] <-  as.character(unique(test_obj$Sample))
  seurat_obj_list[[i]] <- readRDS(obj_name)
}



scrna_GEX <- merge(x=seurat_obj_list[[1]],y=seurat_obj_list[2:length(sample_names)],add.cell.ids = sample_names,project=project_name)

new_cellnames = paste0(project_name,'_',Cells(scrna_GEX))
seurat_obj <- RenameCells(object = scrna_GEX, new.names = new_cellnames)

seurat_obj <- renormalize_recluster(scrna_GEX=seurat_obj,date=date,organism=organism)

saveRDS(seurat_obj, file = paste0('Cycling.SCT.PCA.UMAP.TSNE.CLUST.',project_name,".merged_multisample.",date,".rds"))
