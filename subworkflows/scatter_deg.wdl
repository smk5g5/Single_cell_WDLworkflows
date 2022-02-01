version 1.0

import "../tasks/LinearChain_htmapdegs.wdl" as htmap


workflow scatter_htmap{

	input {
	String docker_image_htmap
	String docker_image_aggr
    String queue_name
    Int mem_gb
    File htmap_rscript
    File aggregate_script
	File inputSamplesFile
	Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)
	}

	scatter (sample in inputSamples) {
		call htmap.aggregate_expression_seurat4_0_3 as aggr {
			input:
			    docker_image=docker_image_aggr,
    			queue_name=queue_name,
    			mem_gb=mem_gb,
    			aggregate_script=aggregate_script,
    			Seurat_rds=sample[2],
    			DEG_file=sample[1],
    			clustering=sample[3],
		}
		call htmap.make_heatmap_DEGs as mkhtmp {
			input:
			    docker_image=docker_image_htmap,
    			queue_name=queue_name,
    			mem_gb=mem_gb,
    			htmap_rscript=htmap_rscript,
    			output_prefix = sample[0],
    			Seurat_rds=sample[2],
    			DEG_file=sample[1],
    			clustering=sample[3],
    			aggregate_exp_rds=aggr.aggregate_rds
		}
	}
}