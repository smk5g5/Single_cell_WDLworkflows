version 1.0

import "../tasks/Doublet_calling.wdl" as doublet_calling
import "../tasks/single_sample_seurat.wdl" as single_filter
import "../tasks/add_doubletinfo.wdl" as add_doubinfo


workflow filter_n_doublets{

	input {
	String Sample_name
	}

	call single_filter.run_seurat_singlesample as seurat_singlesample { 
		input:
		Sample_name=Sample_name,
	}

	call doublet_calling.run_doublet_collection as dc {
		input:
		Sample_name=Sample_name
	}

	call add_doubinfo.run_doublet_collection as add_in {
		input:
		Sample_name=Sample_name
		doublet_file=dc.doublet_results
		input_rds_file=seurat_singlesample.intermed_rds
		output_rds_file=sub(seurat_singlesample.intermed_rds, "\\.rds", ".doublets.rds")
	}

	output {
	    File clust_rds=add_in.seurat_doublet_rds
	}
}