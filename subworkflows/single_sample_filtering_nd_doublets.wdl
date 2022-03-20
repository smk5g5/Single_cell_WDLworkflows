version 1.0

import "../tasks/Doublet_calling.wdl" as doublet_calling
import "../tasks/single_sample_seurat.wdl" as single_filter
import "../tasks/add_doubletinfo.wdl" as add_doubinfo


workflow filter_n_doublets{

	input {
	String Sample_name
	String cellranger_outs_directory
	String docker_image
	String queue_name
	Int mem_gb
	}

	call single_filter.run_seurat_singlesample as seurat_singlesample { 
		input:
		Sample_name=Sample_name,
		cellranger_outs_directory=cellranger_outs_directory,
		docker_image=docker_image,
		mem_gb=mem_gb,
		queue_name=queue_name
	}

	call doublet_calling.run_doublet_collection as dc {
		input:
		Sample_name=Sample_name,
		cellranger_outs_directory=cellranger_outs_directory,
		docker_image=docker_image,
		mem_gb=mem_gb,
		queue_name=queue_name
	}

	call add_doubinfo.add_doublets_metadata as add_in {
		input:
		doublet_file=dc.doublet_results,
		docker_image=docker_image,
		mem_gb=mem_gb,
		queue_name=queue_name,
		input_rds_file=seurat_singlesample.intermed_rds,
	}

	output {
	    File doubmeta_rds=add_in.seurat_doublet_rds
	}
}