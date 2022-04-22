.libPaths( c("/storage1/fs1/allegra.petti/Active/R_libs_scratch/RLibs_4.0.3",.libPaths()) )
library(Seurat)
library(ggplot2)
library(cowplot)
library(dplyr)
library(Matrix)
library(RColorBrewer)
library(ggthemes)
library(GetoptLong)


GetoptLong(
    "seurat_rds=s", "Number of items.",
    "metacol=s", "Name of metadata column to use as Ident.",
    "char_idents=s@", "character vector of idents to remove or keep if inverse is false",
    "inverse=s", "Boolean TRUE or FALSE",
    "append_word=s", "output suffix of rds filtering",
    "organism=s", "Organism name mouse/human.",
    "verbose",  "Print message."
)

# spec = "
# This is an example of using this R script for complex subsetting.

# Usage: Rscript Complex_subset_renorm_recluster.R [options]

# Options:
#   <seurat_rds=s> Path of seurat RDS
#   <metacol=n> Name of metadata column to use as Ident.
#   <idents_vec=v> vector of idents to remove           
#   <inverse=i> Boolean TRUE or FALSE.
#   <output_suffix=s> output suffix of rds filtering.
#   <organism=o> Organism name mouse/human.
# "

# GetoptLong(spec, template_control = list(opt_width = 21))


print("seurat_rds")
print(seurat_rds)
print("metacol")
print(metacol)
print("idents_vec")
print(vec_idents)
print("inverse")
print(inverse)
print("output_suffix")
print(append_word)
print("organism")
print(organism)

# args <- commandArgs(trailingOnly = TRUE)
# if(length(args) < 6) {
#   args <- c("--help")
# }

# seurat_loc <- as.character(args[1])
# sub_col <- as.character(args[2])
# ident_names <- as.character(args[3])
# inverse <- as.character(args[4])
# output_suffix <- as.character(args[5])
# organism <- as.character(args[6])
# # output.stats <- as.character(args[4])
# # output_meta <- as.character(args[5])

# seurat_obj <- readRDS(seurat_loc)

# print(names(seurat_obj@meta.data))

# date = gsub("2022-","22",Sys.Date(),perl=TRUE);
# date = gsub("-","",date);

