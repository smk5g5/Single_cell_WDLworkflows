version 1.0

import "../tasks/clustering_n_pca_simple.wdl" as clus
import "../tasks/single_sample_seurat.wdl" as single_filter


workflow filter_n_cluster{

	input {
	String Sample_name
	}

	call single_filter.run_seurat_singlesample as seurat_singlesample { 
		input:
		Sample_name=Sample_name,
	}

	call clus.run_clustering_n_pca_simple as clus_n_pca {
		input:
		Sample_name=Sample_name,
		String rds_file_path=select_first(seurat_singlesample.intermed_rds)
	}

	output {
	    File clust_rds = clus_n_pca.intermed_rds
	}
}