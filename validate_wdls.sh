java -jar ../womtool-53.1.jar validate ./tasks/scsorter_run.wdl > ../logs/scsorter_run.log
java -jar ../womtool-53.1.jar validate ./tasks/multisample_seurat.wdl > ../logs/multisample_seurat.log
java -jar ../womtool-53.1.jar validate ./tasks/run_cellid.wdl > ../logs/run_cellid.log
java -jar ../womtool-53.1.jar validate ./tasks/linear_chain_adddoublets_seurat.wdl > ../logs/linear_chain_adddoublets_seurat.log
java -jar ../womtool-53.1.jar validate ./tasks/singleR_immune_nonimmune_pipeline.wdl > ../logs/singleR_immune_nonimmune_pipeline.log
java -jar ../womtool-53.1.jar validate ./tasks/plot_qc.wdl > ../logs/plot_qc.log
java -jar ../womtool-53.1.jar validate ./tasks/SingleR_singleref.wdl > ../logs/SingleR_singleref.log
java -jar ../womtool-53.1.jar validate ./tasks/add_doubletinfo.wdl > ../logs/add_doubletinfo.log
java -jar ../womtool-53.1.jar validate ./tasks/clustering_n_pca_simple.wdl > ../logs/clustering_n_pca_simple.log
java -jar ../womtool-53.1.jar validate ./tasks/rpca_integration.wdl > ../logs/rpca_integration.log
java -jar ../womtool-53.1.jar validate ./tasks/htmap_degs.wdl > ../logs/htmap_degs.log
java -jar ../womtool-53.1.jar validate ./tasks/add_cellid_results_to_seurat.wdl > ../logs/add_cellid_results_to_seurat.log
java -jar ../womtool-53.1.jar validate ./tasks/Doublet_calling.wdl > ../logs/Doublet_calling.log
java -jar ../womtool-53.1.jar validate ./tasks/make_gene_featplots.wdl > ../logs/make_gene_featplots.log
java -jar ../womtool-53.1.jar validate ./tasks/LinearChain_htmapdegs.wdl > ../logs/LinearChain_htmapdegs.log
java -jar ../womtool-53.1.jar validate ./tasks/single_sample_seurat.wdl > ../logs/single_sample_seurat.log
java -jar ../womtool-53.1.jar validate ./tasks/subset_renormalize_recluster.wdl > ../logs/subset_renormalize_recluster.log
java -jar ../womtool-53.1.jar validate ./subworkflows/recluster_renorm_rerun_singleR.wdl > ../logs/recluster_renorm_rerun_singleR.log
java -jar ../womtool-53.1.jar validate ./subworkflows/scatter_gather_singleR.wdl > ../logs/scatter_gather_singleR.log
java -jar ../womtool-53.1.jar validate ./subworkflows/single_sample_filtering_nd_doublets.wdl > ../logs/single_sample_filtering_nd_doublets.log
java -jar ../womtool-53.1.jar validate ./subworkflows/scatter_deg.wdl > ../logs/scatter_deg.log
java -jar ../womtool-53.1.jar validate ./subworkflows/scatter_gather_Doubletcall_multisample.wdl > ../logs/scatter_gather_Doubletcall_multisample.log
java -jar ../womtool-53.1.jar validate ./subworkflows/seurat_counts_to_anndata.wdl > ../logs/seurat_counts_to_anndata.log
java -jar ../womtool-53.1.jar validate ./subworkflows/scatter_gather_scsorter.wdl > ../logs/scatter_gather_scsorter.log
java -jar ../womtool-53.1.jar validate ./subworkflows/single_sample_filtering_n_clustering.wdl > ../logs/single_sample_filtering_n_clustering.log
java -jar ../womtool-53.1.jar validate ./subworkflows/scatter_gather_doublet_singlesample.wdl > ../logs/scatter_gather_doublet_singlesample.log
java -jar ../womtool-53.1.jar validate ./subworkflows/scatter_gather_singleR_bysample.wdl > ../logs/scatter_gather_singleR_bysample.log
java -jar ../womtool-53.1.jar validate ./subworkflows/scatter_rpca.wdl > ../logs/scatter_rpca.log
java -jar ../womtool-53.1.jar validate ./pipelines/multisample_seurat_doublet_singler.wdl > ../logs/multisample_seurat_doublet_singler.log
java -jar ../womtool-53.1.jar validate ./pipelines/end_to_end_seurat_singlesample.wdl > ../logs/end_to_end_seurat_singlesample.log
java -jar ../womtool-53.1.jar validate ./pipelines/end_to_end_multisample.wdl > ../logs/end_to_end_multisample.log
java -jar ../womtool-53.1.jar validate ./pipelines/seurat_single_sample_nomerging.wdl > ../logs/seurat_single_sample_nomerging.log
