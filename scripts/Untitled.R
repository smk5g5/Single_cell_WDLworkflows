library('GEOquery')
library(Biobase)
library(oligoClasses)
#Annotation and data import packages
library(ArrayExpress)
library(pd.hugene.1.0.st.v1)
library(hugene10sttranscriptcluster.db)
#Quality control and pre-processing packages
library(oligo)
library(arrayQualityMetrics)
#Analysis and statistics packages
library(limma)
library(topGO)
library(ReactomePA)
library(clusterProfiler)
library(pheatmap)
library(hta20transcriptcluster.db)
library(limma)

library("RColorBrewer")
library(tidyverse)
library(ComplexHeatmap)
library(circlize)
library(genefilter)

#This script will be updated to work with WDL

makecomplexheatmap_bysymbols2 <- function(aggr_res,sel_genes,annot_df,des_col,gsmid,no_of_genes='all',myhite=50,rowfont,htmfont,mywid=30) {
  set.seed(100)
  gse_df <- aggr_res$combined_df
  sub_gsedf <- gse_df[gse_df$Gene_Symbol %in% sel_genes,]
  eset <- aggr_res$combined_eset
  
  qual_cols = brewer.pal.info[brewer.pal.info$category == "qual", ]
  qual_cols <- qual_cols[qual_cols$colorblind==T,]
  qual_cols
  col_vector = unlist(mapply(brewer.pal, qual_cols$maxcolors, rownames(qual_cols)))
  
  eset_ligs <- eset[rownames(eset) %in% sub_gsedf$Gene_Symbol,]
  
  gene_len <- length(rownames(eset_ligs))
  base_mean = rowMeans(eset_ligs)
  mat_scaled = t(scale(t(eset_ligs)))
  
  index <- match(colnames(eset_ligs),merge_annot_final$gsmid)
  groupdf <- merge_annot_final[['Type']][index]
  subtype_df <- merge_annot_final[['sub_type']][index]
  #hearing_df <- merge_annot_final[['Hearing.loss']][index]
  #nf2_df <-  merge_annot_final[['nf2']][index]
  col_cols <- toupper(c('#f1a340','#998ec3'))
  #names(col_cols) <- c('Control','Tumor')
  names(col_cols) <- c('Vestibular nerve tissue','Vestibular schwannoma tumor')
  
  subtypes_cols <- sample(col_vector,length(unique(merge_annot_final$sub_type)))
  names(subtypes_cols) <- unique(merge_annot_final$sub_type)
  
  #hearing_loss_cols <- c('#ef304c','#14c99c','#EEEEEE')
  #names(hearing_loss_cols) <- c('YES','NO','NA')
  
  #nf2_cols <-  c('#ab002b','#00856F','#4D4D4D')
  #names(nf2_cols) <- c('YES','NO','NA')
  #Hearing_loss=hearing_loss_cols,NF2=nf2_cols  
  #Hearing_loss=hearing_df,NF2=nf2_df
  ha = HeatmapAnnotation(df = data.frame(type = groupdf,subtype=subtype_df),col = list(type = col_cols,subtype=subtypes_cols),gp=gpar(fontsize = 4))
  
  #des_col = 'Control - Tumor'
  sub_gsedf$genereg[sub_gsedf[[des_col]]==-1] <- 'Downregulation'
  sub_gsedf$genereg[sub_gsedf[[des_col]]==1] <- 'Upregulation'
  sub_gsedf$genereg[sub_gsedf[[des_col]]==0] <- 'No_change'
  
  index <- match(rownames(mat_scaled),sub_gsedf$Gene_Symbol)
  # rownames(mat_scaled) <- sub_gsedf$Gene_Symbol[index]
  selvals <- sub_gsedf$genereg[index]
  
  var_genes <- sub_gsedf$variable_Genes[index]
  
  row_ha = rowAnnotation(diffexp=selvals,variable_gene=var_genes,col = list(diffexp = c("Downregulation" = "red", "Upregulation" = "green", "No_change" = "grey"),variable_gene=c("YES"="#dd2c00","NO"="#8fe7ff")),width = unit(1, "cm"))
  ht <- Heatmap(mat_scaled , name=sprintf("Z-score normalized gene expression matrix %s",gsmid),
                col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
                top_annotation = ha,right_annotation = row_ha,
                row_names_gp = grid::gpar(fontsize = rowfont, fontface = "bold"),
                show_row_names = T, show_column_names = F,
                show_row_dend = F,show_column_dend = F)
  heatmap_plot = draw(ht,legend_title_gp = gpar(fontsize = htmfont, fontface = "bold"),legend_grid_width = unit(1, "cm"), legend_grid_height = unit(1, "cm"))
  png(filename = sprintf("heatmap_%s_%s.png",gsmid,no_of_genes),
      width = mywid, height =myhite , units = "in",
      bg = "white",res=300)
  # Heatmap(mat, name = "mat", cluster_rows = row_dend)
  print(heatmap_plot)
  dev.off()
  myrow_ord <- row_order(ht)
  # print(head(myrow_ord))
  # print(length(myrow_ord))
  # print(length(rownames(mat_scaled)))
  #myrow_dend <- row_dend(ht)
  # mycol_ord <- column_order(heatmap_plot)
  # mycol_dend <- column_dend(heatmap_plot)
  # # row_order = sort(rownames(mat)), 
  # # column_order = sort(colnames(mat))
  # # row_km = 2)
  # heatmap_plot <- Heatmap(mat_scaled , name=sprintf("Z-score normalized gene expression matrix %s",gsmid),
  #                         col=colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
  #                         top_annotation = ha,right_annotation = row_ha,
  #                         row_names_gp = grid::gpar(fontsize = 10),
  #                         show_row_names = T, show_column_names = F,
  #                         show_row_dend = T,show_column_dend = T,row_order = unlist(myrow_ord),
  #                         column_order = mycol_ord,
  #                         row_dend_width = unit(80, "mm"),
  #                         column_dend_height = unit(80, "mm"),row_km =20,row_km_repeats=10)
  # png(filename = sprintf("heatmap_ligands_row_col_clus_%s_%s.png",gsmid,no_of_genes),
  #     width = 30, height = 35, units = "in",
  #     bg = "white",res=300)
  # print(heatmap_plot)
  # dev.off()
  ordered_gene_list <- rownames(mat_scaled)[unlist(myrow_ord)]
  return(list(ht=ht,scl_mat=mat_scaled,ord_row=myrow_ord,col_ord=column_order(ht),ord_gene_list=ordered_gene_list))
}