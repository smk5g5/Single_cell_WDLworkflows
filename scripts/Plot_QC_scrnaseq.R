plot_qc_metrics <- function(seurat_obj,filename,mt_pat="^MT-") {
seurat_obj$log10GenesPerUMI <- log10(seurat_obj$nFeature_RNA) / log10(seurat_obj$nCount_RNA)
seurat_obj$mitoRatio <- PercentageFeatureSet(object = seurat_obj, pattern = mt_pat)
seurat_obj$mitoRatio <- seurat_obj@meta.data$mitoRatio / 100
metadata <- seurat_obj@meta.data
metadata$cells <- rownames(metadata)

metadata <- metadata %>%
  dplyr::rename(seq_folder = orig.ident,
                nUMI = nCount_RNA,
                nGene = nFeature_RNA)

Ncells <- metadata %>% 
  ggplot(aes(x=Sample, fill=Sample)) + 
  geom_bar() +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells")

cell_dens <- metadata %>% 
  ggplot(aes(color=Sample, x=nUMI, fill= Sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  ylab("Cell density") +
  geom_vline(xintercept = 500)+ggtitle("Cell density")

genes_per_cell <- metadata %>% 
  ggplot(aes(color=Sample, x=nGene, fill= Sample)) + 
  geom_density(alpha = 0.2) + 
  theme_classic() +
  scale_x_log10() + 
  geom_vline(xintercept = 300) + ggtitle("Genes per cell")

box_genespercell <- metadata %>% 
  ggplot(aes(x=Sample, y=log10(nGene), fill=Sample)) + 
  geom_boxplot() + 
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
  theme(plot.title = element_text(hjust=0.5, face="bold")) +
  ggtitle("NCells vs NGenes")

umibygenesmito <- metadata %>% 
  ggplot(aes(x=nUMI, y=nGene, color=mitoRatio)) + 
  geom_point() + 
  scale_colour_gradient(low = "gray90", high = "black") +
  stat_smooth(method=lm) +
  scale_x_log10() + 
  scale_y_log10() + 
  theme_classic() +
  geom_vline(xintercept = 500) +
  geom_hline(yintercept = 250) +
  facet_wrap(~Sample) + ggtitle("UMI vs Genes colored by mito ratio")


mitoratio_dens <- metadata %>% 
  ggplot(aes(color=Sample, x=mitoRatio, fill=Sample)) + 
  geom_density(alpha = 0.2) + 
  scale_x_log10() + 
  theme_classic() +
  geom_vline(xintercept = 0.2) + ggtitle("mitoratio density")

genes_per_umi <- metadata %>%
  ggplot(aes(x=log10GenesPerUMI, color = Sample, fill=Sample)) +
  geom_density(alpha = 0.2) +
  theme_classic() +
  geom_vline(xintercept = 0.8) + ggtitle("Genes per UMI")

pdf(filename)
print(Ncells)
print(cell_dens)
print(genes_per_cell)
print(box_genespercell)
print(umibygenesmito)
print(mitoratio_dens)
print(genes_per_umi)
dev.off()
}

plot_all_samples <- function(seurat_obj){
  all_samples <- unique(seurat_obj$Sample)
  plots=list();
  color_codes <- brewer.pal(length(all_samples),"Spectral")
  names(color_codes) <- all_samples
  jpeg(sprintf("%s/UMAP.50.ColorbySample%s.%s.jpg", output.stats, control, date), width = 10, height = 8, units="in", res=300);
  p2 <- DimPlot(object = scrna, reduction = "umap", group.by = "Sample", pt.size=0.1,cols = color_codes)
  print(p2)
  dev.off()
  for (i in 1:length(all_samples)) {
    samp <- all_samples[i]
    subcells <- Cells(seurat_obj)[(which(seurat_obj[['Sample']][['Sample']] == samp))];
    subcolors <- color_codes[samp]
    plots[[i]] <- DimPlot(object = seurat_obj, reduction = "umap", group.by = "Sample", cells.highlight=subcells, cols.highlight=subcolors, sizes.highlight = 1, pt.size=0.9) + theme(axis.title.x=element_blank(),axis.title.y=element_blank(),axis.text.x=element_blank(),axis.text.y=element_blank(),axis.ticks.x=element_blank(),axis.ticks.y=element_blank());
  }
  jpeg(sprintf("%s/Combinedplots_ColorbySample%s.%s.jpg",output.stats, control, date), width = 50, height = 50, units="in", res=300);
  print(CombinePlots(plots, ncol=4, legend="none")); # or 4 for 12,16,20
  dev.off()
  }

