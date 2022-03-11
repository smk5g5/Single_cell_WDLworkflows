version 1.0

import "../tasks/Doublet_calling.wdl" as doublet

workflow scatter_doublet{

  input {
    String docker_image
    String queue_name
    Int mem_gb
    File doublet_calling_script
    File inputSamplesFile
    Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)
  }

  scatter (sample in inputSamples) {
    call doublet.run_doublet_collection as run_doublet {
      input:
          docker_image=docker_image,
          queue_name=queue_name,
          mem_gb=mem_gb,
          doublet_calling_script=doublet_calling_script,
          cellranger_outs_directory=sample[1],
          Sample_name=sample[0]
  }

}
}