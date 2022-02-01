library(Seurat)
library(ggplot2)
library(stringr)
library(ComplexHeatmap)
library(circlize)
library(tidyverse)
library(dplyr)
library(clusterProfiler)
library(org.Hs.eg.db)
library(RColorBrewer)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 3) {
  args <- c("--help")
}

Seurat_file <- as.character(args[1])
DEG_file <- as.character(args[2])
prefix <- as.character(args[3])
clustering <- as.character(args[4])

seurat_object <- readRDS(Seurat_file)

DEGs <- read.table(DEG_file,header=T)

date = gsub("2021-","21",Sys.Date(),perl=TRUE);
date = gsub("-","",date);


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


make_heatmap_compdo <- function(seurat_object,clustering,sel_genes,anno_genes,date,prefix){

set.seed(100)

clust.means.norm <- get_avg_scaledexp(seurat_object = seurat_object,clustering=clustering,sel_genes=sel_genes)

colors = unique(c(seq(-1,-0.1,length=27),seq(-0.1,0.1,length=24),seq(0.1,1,length=27)));
my.palette = colorRampPalette(rev(brewer.pal(n=7, name="RdYlBu")))(n=75);

myhtmp_phtmp <- ComplexHeatmap::pheatmap(clust.means.norm,cluster_rows = T,cluster_cols = T,show_rownames = T,column_title = "Average expression by cluster",show_colnames=T,treeheight_row=0,treeheight_col=0,
                                          border_color='black', scale="none",breaks=colors,color=my.palette,fontsize=30,fontsize_row=16,fontsize_col=30,
                                          cellwidth = ncol(clust.means.norm)*unit(10, "mm"), cellheight = nrow(clust.means.norm)*unit(5, "mm"))

fig_ht <- min(125,nrow(clust.means.norm)*unit(5, "mm"))
fig_wd <- min(160,2 * ncol(clust.means.norm)*unit(5, "mm"))


png(filename=sprintf("%s_RNAassay.%s.aggregatedhtmap.%s.png",prefix,clustering,date), units="cm", res=300, height=fig_ht, width=fig_wd);
print(myhtmp_phtmp)
dev.off()

gene_indices <- match(anno_genes,rownames(clust.means.norm))

# row_ha = rowAnnotation(diffexp=selvals,col = list(diffexp = c("Downregulation" = "red", "Upregulation" = "green", "No_change" = "grey")),width = unit(2, "cm"))
gene_annots = rowAnnotation(selgenes = anno_mark(at = gene_indices, labels = anno_genes))


myhtmp_phtmp_anno <- ComplexHeatmap::pheatmap(clust.means.norm,cluster_rows = T,cluster_cols = T,show_rownames = F,right_annotation = gene_annots,column_title = "Average expression by cluster",show_colnames=T,treeheight_row=0,treeheight_col=0,
                                          border_color='black', scale="none",breaks=colors,color=my.palette,fontsize=30,fontsize_row=16,fontsize_col=30,
                                          cellwidth = ncol(clust.means.norm)*unit(10, "mm"), cellheight = nrow(clust.means.norm)*unit(5, "mm"))

heatmap_plot = draw(myhtmp_phtmp_anno,legend_title_gp = gpar(fontsize = 10, fontface = "bold"),legend_grid_width = unit(1, "cm"), legend_grid_height = unit(1, "cm"))
png(filename = sprintf("%s_RNAassay.%s.aggregatedhtmap_top5DEGpercluster_annot.%s.png",prefix,clustering,date),
    width = 22, height =15 , units = "in",
    bg = "white",res=600)
# Heatmap(mat, name = "mat", cluster_rows = row_dend)
print(heatmap_plot)
dev.off()


doheatmap_roworder_genes <- rownames(clust.means.norm)[unname(unlist(row_order(heatmap_plot)))]
doheatmap_col_ordr <- colnames(myscaled_htmapdata)[column_order(heatmap_plot)]


levels(seurat_object) <- doheatmap_col_ordr
mini <- subset(seurat_object, features=doheatmap_roworder_genes);
levels(mini) <- doheatmap_col_ordr


png(filename=sprintf("%s_DoHeatmap.By%s.RNAassay.%s.png",prefix,clustering,date), width=35, height=50, units="in", res=300);
hm <- DoHeatmap(mini, features=doheatmap_roworder_genes, slot="data", disp.min=-1.5, disp.max=7, group.by="ident", group.bar=TRUE) + scale_fill_gradientn(colors = c("blue", "white", "red"))
print(hm);
dev.off();

saveRDS(list(phtmap=myhtmp_phtmp,htmap_draw=heatmap_plot,roword=doheatmap_roworder_genes,col_ord=doheatmap_col_ordr,doheatmap=hm),file=sprintf("%s_RNAassay.%s.aggregatedhtmap.%s.rds",prefix,clustering,date))
}


# sel_genes <- c("S100B","NCAM1","PMP22","CDH19","ITGB8","SLIT2","ADAM23","CADM2","NCAM2","NRXN1","SCN7A","NRXN3")
# gene_indices <- match(sel_genes,rownames(mat_scaled))
# gene_annots = rowAnnotation(selgenes = anno_mark(at = gene_indices, labels = sel_genes))

# gene_annots = rowAnnotation(selgenes = anno_mark(at = gene_indices, labels = sel_genes))

# ht <- Heatmap(mat_scaled , name=sprintf("Z-score normalized gene expression matrix %s",'GSE141801'),
#               col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
#               top_annotation = ha,left_annotation = row_ha,right_annotation = gene_annots,
#               row_names_gp = grid::gpar(fontsize = 10, fontface = "bold"),
#               show_row_names = F, show_column_names = F,
#               show_row_dend = F,show_column_dend = F,rect_gp = gpar(col = "black", lwd = 1),width = ncol(mat_scaled)*unit(5, "mm"), 
#               height = nrow(mat_scaled)*unit(5, "mm"))
# heatmap_plot = draw(ht,legend_title_gp = gpar(fontsize = 10, fontface = "bold"),legend_grid_width = unit(1, "cm"), legend_grid_height = unit(1, "cm"))
# png(filename = sprintf("heatmap_%s_%s.png",'GSE141801','top50_DEGs_nmSC-like'),
#     width = 22, height =10 , units = "in",
#     bg = "white",res=600)
# # Heatmap(mat, name = "mat", cluster_rows = row_dend)
# print(heatmap_plot)
# dev.off()



DEGs_sig <- DEGs[DEGs$p_val_adj <= 0.05,]

top25_DEGs_sig <- DEGs_sig %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC))%>% 
  dplyr::slice(1:25)

top5_DEGs_sig <- DEGs_sig %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC))%>% 
  dplyr::slice(1:5)

make_heatmap_compdo(seurat_object = seurat_object,clustering=clustering,sel_genes=unique(top25_DEGs_sig$gene),anno_genes=unique(top5_DEGs_sig$gene),date=date,prefix=prefix)

DEG_marks_list <- list()
for(i in unique(top25_DEGs_sig$cluster)){
  subclus_deg <- subset(top25_DEGs_sig,cluster==i)
  DEG_marks_list[[paste0('Cluster_',i)]] <- subclus_deg$gene
}

compcluster_out <- compareCluster(geneCluster = DEG_marks_list,OrgDb = org.Hs.eg.db, keyType="SYMBOL", fun='enrichGO',ont="BP")

jpeg(sprintf("%s.enrichment_summary_BP_top25DEGs_%s_showtopCategory.%s.jpg",prefix,cluster,date), width = 15, height = 20, units="cm", res=600);
dotplot(compcluster_out) + theme(axis.text.x = element_text(angle = 90))
dev.off()


jpeg(sprintf("%s.enrichment_summary_BP_top25DEGs_%s_showCategory50.%s.jpg",prefix,cluster,date), width = 20, height = 50, units="in", res=300);
dotplot(compcluster_out,showCategory=50) + theme(axis.text.x = element_text(angle = 90))
dev.off()


cell_marker_data <- vroom::vroom('http://bio-bigdata.hrbmu.edu.cn/CellMarker/download/Human_cell_markers.txt')

## instead of `cellName`, users can use other features (e.g. `cancerType`)
cells <- cell_marker_data %>%
    dplyr::select(cellName, geneSymbol) %>%
    dplyr::mutate(geneSymbol = strsplit(geneSymbol, ', ')) %>%
    tidyr::unnest(cols = c(geneSymbol))

jpeg(sprintf("%s.celltype_enrichment_top25DEGs_%s.%s.jpg",prefix,cluster,date), width = 15, height = 20, units="cm", res=600);
dotplot(celltype_enrichment_merged) + theme(axis.text.x = element_text(angle = 90))
dev.off()

celltype_enrichment_list <-  lapply(DEG_marks_list,enricher,TERM2GENE = cells)
saveRDS(list(celltype_enrichment=celltype_enrichment,BP_enrichment=compcluster_out),file=sprintf("%s_RNAassay.%s.enrichment_compclus_celltype.%s.rds",prefix,clustering,date))


jpeg("celltype_enrichment_top25DEGs_schwann_scrnaonly.jpg", width = 25, height = 30, units="cm", res=600);
dotplot(scrna_celltype_enrichment_merged) + theme(axis.text.x = element_text(angle = 90))
dev.off()

