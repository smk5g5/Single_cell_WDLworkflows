version 1.0

import "../tasks/htmap_degs.wdl" as htmap


workflow scatter_htmap{

	input {
	String docker_image
    String queue_name
    Int mem_gb
    File htmap_rscript
	File inputSamplesFile
	Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)
	}

	scatter (sample in inputSamples) {
		call htmap.make_heatmap_DEGs as mk_htmap {
			input:
			    docker_image=docker_image,
    			queue_name=queue_name,
    			mem_gb=mem_gb,
    			htmap_rscript=htmap_rscript,
    			Seurat_rds=sample[1],
    			output_prefix=sample[0],
    			DEG_file=sample[2]
    			clustering=sample[2]
	}

}
}