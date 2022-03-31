library(CelliD)
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)
if(length(args) < 3) {
  args <- c("--help")
}

Seurat_file <- as.character(args[1])
marker_gene_list_rds <- as.character(args[2])
marker_file_reference_name <- as.character(args[3])

seurat_object <- readRDS(Seurat_file)
marker_list <- readRDS(marker_gene_list_rds)

date = gsub("2022-","22",Sys.Date(),perl=TRUE);
date = gsub("-","",date);

#cellid docker image we are using is Seurat_4.1.0
#we are basically taking seurat 4.0.3 object and making cellid predictions using Seurat_4.1.0
# seurat_object <- readRDS('/storage1/fs1/allegra.petti/Active/Users/khan.saad/WDL_pipelines/seurat_counts_to_anndata/0d1c2574-28a1-49e1-88ce-537ad2e57312/call-Get_seurat_counts/inputs/-1835829567/Cycling.SCT.PCA.UMAP.TSNE.CLUST.211216.rds')

#only use protein coding genes
seurat_object_prot <- seurat_object[rownames(seurat_object) %in% HgProteinCodingGenes,]

#To perform MCA dimensionality reduction, the command RunMCA is used:

seurat_object_prot <- RunMCA(seurat_object_prot)

# marker_list <- readRDS('/scratch1/fs1/allegra.petti/khan.saad/WDL_workflow/Single_cell_WDLworkflows/essential_inputs/Wang_CancerCell_list.rds')

# DimPlotMC(Baron, reduction = "mca", group.by = "cell.type", features = c("CTRL", "INS", "MYZAP", "CDH11"), as.text = T) + ggtitle("MCA with some key gene markers")

######################################################
#GBM specific celltype gene lists are in essential inputs rds files
#for other celltype marker lists panglao can be used when necessary
######################################################
# download all cell-type gene signatures from panglaoDB
# panglao <- read_tsv("https://panglaodb.se/markers/PanglaoDB_markers_27_Mar_2020.tsv.gz")

# # restricting the analysis to pancreas specific gene signatues
# panglao_pancreas <- panglao %>% filter(organ == "Pancreas")

# # restricting to human specific genes
# panglao_pancreas <- panglao_pancreas %>%  filter(str_detect(species,"Hs"))

# # converting dataframes into a list of vectors, which is the format needed as input for CellID
# panglao_pancreas <- panglao_pancreas %>%  
#   group_by(`cell type`) %>%  
#   summarise(geneset = list(`official gene symbol`))
# pancreas_gs <- setNames(panglao_pancreas$geneset, panglao_pancreas$`cell type`)


# Performing per-cell hypergeometric tests against the gene signature collection

# n.features: integer of top n features to consider for hypergeometric
 # minSize: minimum number of overlapping genes in geneset and

#there are a number of parameters that can be changed with this i.e no of top features to consider for hypergeometric test
#minimum size of overlap to consider for calculating hypergeometric test. The default is 10 but I left it to 5
#here as 10 gives error with most gene lists. e.g.
#   |++++++++++++++++++++++++++++++++++++++++++++++++++| 100% elapsed=12s
# Error in RunCellHGT.Seurat(seurat_object_prot, pathways = marker_list,  :
#   All pathways have less than 10 features in common with the data
#if we use top 200 features only lot of the cells are assigned unassigned so using top 2000 features instead

#using 10 as the minsize or least number of genes in a marker list
#as some lists in glioblastoma/brain celltypes in essential inputs have only 5 genes.
min_size <- min(10,min(unlist(lapply(marker_list,length))))

HGT_gs <- RunCellHGT(seurat_object_prot, pathways = marker_list, dims = 1:50, n.features = 2000,minSize=min_size)

# For each cell, assess the signature with the lowest corrected p-value (max -log10 corrected p-value)
cellid_gs_prediction <- rownames(HGT_gs)[apply(HGT_gs, 2, which.max)]

# For each cell, evaluate if the lowest p-value is significant
cellid_gs_prediction_signif <- ifelse(apply(HGT_gs, 2, max)>2, yes = cellid_gs_prediction, "unassigned")

#save cellid results for adding them to seurat object assay and metadata.

saveRDS(cellid_gs_prediction_signif, file = paste0('CellID_final_predictions.',marker_file_reference_name,".",date,".rds"))
saveRDS(HGT_gs, file = paste0('CellID_pval_matrix.',marker_file_reference_name,".",date,".rds"))


