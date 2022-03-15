version 1.0

import "../tasks/SingleR_singleref.wdl" as singleR


workflow scatter_singleR{

	input {
	String docker_image
    String queue_name
    Int mem_gb
    File singleR_singleref_rscript
    File seurat_rds
	File inputSamplesFile
	Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)
	}

	scatter (sample in inputSamples) {
		call singleR.run_singleR_singleref as run_singleR {
			input:
			    docker_image=docker_image,
    			queue_name=queue_name,
    			mem_gb=mem_gb,
    			singleR_singleref_rscript=singleR_singleref_rscript,
    			Seurat_rds=seurat_rds,
    			reference_name=sample[0],
    			singleR_ref_rds=sample[2],
				label_column_name=sample[1]
	}}
}