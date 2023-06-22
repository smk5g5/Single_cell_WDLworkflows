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
library(RColorBrewer)
library(gridExtra)


args <- commandArgs(trailingOnly = TRUE)
# if(length(args) < 3) {
#   args <- c("--help")
# }

Seurat_file <- as.character(args[1])
singleR_tsv <- as.character(args[2])
#sample_name <- as.character(args[3])
singleR_files <- as.character(args[3:length(args)])

date = gsub("2023-","23",Sys.Date(),perl=TRUE);
date = gsub("-","",date);

Add_singleR_scores_to_seuratassay_singleref <- function(singleR_obj,seurat_obj,reference_name){
  singleR_scores <- singleR_obj[['scores']]
  rownames(singleR_scores) <- rownames(singleR_obj)
  index <- match(Cells(seurat_obj),rownames(singleR_scores))
  singleR_scores <- singleR_scores[index,]
  singleR_score_assay <- CreateAssayObject(data = t(as.matrix(singleR_scores)))
  seurat_obj[[sprintf("%s.singleR.scores", reference_name)]] <- singleR_score_assay
return(seurat_obj)
}

Add_singleR_preds_to_seuratmeta <- function(singleR_obj,seurat_obj,reference_name){
  #make named vector for cells which are not malignant and not immune

  singleR_preds <- singleR_obj$pruned.labels
  names(singleR_preds) <- rownames(singleR_obj)
  singleR_preds[which(is.na(singleR_preds))] = "No.Prediction"
  #assign celltypes to seurat object to metadata name
  index <- match(Cells(seurat_obj),names(singleR_preds))
  seurat_obj[[sprintf("singleR_results_%s",reference_name)]] <- singleR_preds[index]
  return(seurat_obj)  
}

plot_singleRhca <- function(seurat_obj,meta_celltype_name,ref_name,date) {
  Idents(seurat_obj) <- meta_celltype_name
  lineage.found.all = sort(unique(Idents(object=seurat_obj))); # all found lineages
  n = length(lineage.found.all)

  #  lineage.found = lineage.found.all[which(lineage.found.all != "NA")]; # found lineag
  #es excluding no prediction
  # gray.index = which(lineage.found.all == "No.Prediction")
  # color with rainbow colors:
  qual_cols = brewer.pal.info[brewer.pal.info$category == "qual", ]
  qual_cols = qual_cols[qual_cols$colorblind=='TRUE',]
  col_vector = unlist(mapply(brewer.pal, qual_cols$maxcolors, rownames(qual_cols)))
  # cell.colors = rainbow(n, s=0.6, v=0.9);
  if(n<=length(col_vector)){
  cell.colors <- sample(col_vector, n)
  names(cell.colors) <- lineage.found.all
  }else{
  	cell.colors = rainbow(n,s=0.6, v=0.9);
  	names(cell.colors) <- lineage.found.all
  }
  # cell.colors["No.Prediction"] <- 'gray52'
  cell.colors["No.Prediction"] <- 'gray88'
  # cell.colors.temp = rainbow.colors[which(lineage.list %in% lineage.found)]; # ensure that cell
  # names(cell.colors.temp) <- lineage.list[which(lineage.list %in% lineage.found)]
  # #type colors are constant - this does not include gray for no prediction
  # if (length(gray.index) > 0) {
  #   cell.colors = c(cell.colors.temp[1:gray.index-1],"gray",cell.colors.temp[gray.index:length(cell.colors.temp)])
  #   ; # insert gray at "No.Prediction" position
  #   names(cell.colors) <- names(cell.colors.temp)
  # } else {
  #   cell.colors = cell.colors.temp;
  #   names(cell.colors) <- names(cell.colors.temp)
  # }
  # jpeg(sprintf("%s/%s_singleR_celltype_%s.jpg",output.dir,outfile,ref_name), width = 20, height = 14, units="in", res=300);
  p2 <- DimPlot(object = seurat_obj, reduction = "umap", group.by = meta_celltype_name, cols = cell.colors, pt.size=2)+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  ggsave(sprintf("Dimplot_SingleR_%s_%s_%s.png",ref_name,sample_name,date),plot = p2, width = 30, height = 30, units = "in",dpi = 300,device = "png",scale = 1)
  # print(p2);
  # dev.off();
  # jpeg(sprintf("%s/%s_celltype_%s_labeled.jpg",output.dir,outfile,ref_name), width = 20, height = 14, units="in", res=300);
  p2 <- DimPlot(object = seurat_obj, reduction = "umap", group.by = meta_celltype_name, cols = cell.colors, pt.size=2,label=T,label.size = 5)+ theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  ggsave(sprintf("Dimplot_SingleR_lab_%s_%s_%s.png",ref_name,sample_name,date),plot = p2, width = 30, height = 30, units = "in",dpi = 300,device = "png",scale = 1)
  # print(p2);
  # dev.off()
  plots=list();
  for (i in 1:length(lineage.found.all)) {
    samp = lineage.found.all[i];
    subcells <- Cells(seurat_obj)[(which(seurat_obj[[meta_celltype_name]][[meta_celltype_name]] == samp))];
    subcolors <- cell.colors[samp];
    plots[[i]] <- DimPlot(object = seurat_obj, reduction = "umap", group.by = meta_celltype_name, cells.highlight=subcells, cols.highlight=subcolors,sizes.highlight=0.1) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank())+ggtitle(samp)+theme(plot.title = element_text(hjust = 0.5))
  }
  outfile <- sprintf("Dimplots_by_celltype_%s.%s.pdf",ref_name,sample_name,date);
  ml <- marrangeGrob(plots, nrow=2, ncol=2)
  ggsave(outfile,ml,width=10, height=10)
  return(cell.colors)
  }


seurat_object <- readRDS(Seurat_file)

input_df <- read.table(singleR_tsv,sep="\t",header=FALSE)
colnames(input_df) <- c('Reference_name','label_column_name','reference_rds')


single_R_preds <- list()

for(i in 1:nrow(input_df)){
 print(input_df$Reference_name[i])
 print(singleR_files[grep(input_df$Reference_name[i],singleR_files)])
single_R_preds <-  readRDS(singleR_files[grep(input_df$Reference_name[i],singleR_files)])
 }

for(i in names(single_R_preds)){
seurat_object <- Add_singleR_scores_to_seuratassay_singleref(singleR_obj=single_R_preds,seurat_obj=seurat_object,reference_name=Reference_name)
seurat_object <- Add_singleR_preds_to_seuratmeta(singleR_obj=single_R_preds,seurat_obj=seurat_object,reference_name=Reference_name)
plot_singleRhca(seurat_obj=seurat_object,meta_celltype_name=sprintf("singleR_results_%s",Reference_name),ref_name=Reference_name,date=date)
}

output_file <- paste0(gsub('\\.[0-9]*.rds$','',basename(Seurat_file)),".single_R_preds.",date,".rds")

saveRDS(seurat_object, file = output_file)

