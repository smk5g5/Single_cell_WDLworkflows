version 1.0

import "../subworkflows/scatter_gather_singleR.wdl" as scatter_gather_singleR

workflow run_singleR_by_sample{

  input {
  String docker_image
    String queue_name
    Int mem_gb
    File singleR_refsFile
  File inputSamplesFile
  Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)
  }

  scatter (sample in inputSamples) {
    call scatter_gather_singleR.scatter_gather_singleR as scat_gath_singleR {
      input:
      docker_image=docker_image,
      queue_name=queue_name,
      mem_gb=mem_gb,
      seurat_rds=sample[0],
      inputSamplesFile=singleR_refsFile
    }
  }
}