.libPaths( c("/storage1/fs1/allegra.petti/Active/R_libs_scratch/RLibs_4.0.3",.libPaths()) )
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(Matrix)
library("sctransform");
library("dplyr");
library("RColorBrewer");
library("ggthemes");
library("ggplot2");
library("cowplot");
library(tidyverse)

library(ggplot2)
library(gridExtra)
#source('./scripts/Plot_QC_scrnaseq.R')

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 4) {
  args <- c("--help")
}

#
seurat_rds <- as.character(args[1])
gene_lists <- as.character(args[2])
gene_column <- as.character(args[3])
type_column <- as.character(args[4])


control <- 'Cycling'

date = gsub("2022-","22",Sys.Date(),perl=TRUE);
date = gsub("-","",date);

glists.raw <- read.table(gene_lists, sep=",",row.names=NULL,header=TRUE,as.is=TRUE); # gene lists
glists <- glists.raw[which(glists.raw[[gene_column]] %in% rownames(scrna_GEX)),]; # filtered genelists

for(i in unique(glists[[type_column]])) {
  j=gsub(" ","_",i);
  j=gsub("/","_",j);
  sub_glists <- glists[glists[[type_column]]==i,]
  genesToPlot = sub_glists[[gene_column]]
  ng = length(genesToPlot); # number of genes
  outfile = sprintf("umap.%s.%s.%s.pdf",j, control, date);
  if(ng==1){
  fp <- FeaturePlot(object = scrna_GEX, features = genesToPlot, cols = c("gray","red"), ncol=2, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  ggsave(outfile,plot=fp,width = 5, height = 5) + ggtitle(sprintf("%s (%s)",k,j))
  }  else{
    gplotlist=list();
    for(k in genesToPlot){
      gplotlist[[k]] <- FeaturePlot(object = scrna_GEX, features = k, cols = c("gray","red"), ncol=2, reduction = "umap") + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank())
    }
    ml <- marrangeGrob(gplotlist, nrow=2, ncol=2)
    ggsave(outfile,ml,width=10, height=10)
  }
}