version 1.0

import "../tasks/clustering_n_pca_simple.wdl" as clus
import "../tasks/single_sample_seurat.wdl" as single_filter


workflow filter_n_cluster{

	input {
	String cellranger_outs_directory
	String Sample_name
	File seurat_singlesample_rscript
	File clustering_script
	String tirosh_file_path
	Float clustering_res
	Int nPC
	}

	call single_filter.run_seurat_singlesample as seurat_singlesample { 
		input:
		cellranger_outs_directory=filter_n_cluster.cellranger_outs_directory,
		Sample_name=filter_n_cluster.Sample_name,
		seurat_singlesample_rscript=filter_n_cluster.seurat_singlesample_rscript
	}

	call clus.run_clustering_n_pca_simple as clus_n_pca {
		input:
		clustering_script=filter_n_cluster.clustering_script,
		Sample_name=filter_n_cluster.Sample_name,
		rds_file_path=seurat_singlesample.intermed_rds
	}

	output {
	    File clust_rds = clus_n_pca.intermed_rds
	}
}