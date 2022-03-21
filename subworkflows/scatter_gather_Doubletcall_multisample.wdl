version 1.0

import "../tasks/Doublet_calling.wdl" as doublet

workflow scatter_doublet{

  input {
    String docker_image
    String queue_name
    Int mem_gb
    File doublet_calling_script
    File merged_seurat_rds
    File merge_doublet_calls_in_seurat_script
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

 call add_doublets_metadata_tomultisample_seurat {
      input: 
        docker_image=docker_image, 
        queue_name=queue_name,
        mem_gb=mem_gb,
        merge_doublet_calls_in_seurat_script=merge_doublet_calls_in_seurat_script,
        inputSamplesFile=inputSamplesFile,
        multisample_seurat_rds=merged_seurat_rds,
        doublet_files=run_doublet.doublet_results
    }

  output {
  File seurat_doublet_rds = add_doublets_metadata_tomultisample_seurat.seurat_doublet_rds
  }

}

task add_doublets_metadata_tomultisample_seurat{

  input {
    String docker_image
    String queue_name
    Int mem_gb
    File merge_doublet_calls_in_seurat_script
    File inputSamplesFile
    File multisample_seurat_rds
    String split_by
    Array[String] doublet_files
  }

   command <<<
    Rscript ~{merge_doublet_calls_in_seurat_script} ~{multisample_seurat_rds} ~{inputSamplesFile} ~{split_by} ~{sep=" " doublet_files}
    >>>

  runtime {
    docker : docker_image
    memory: mem_gb + " GB"
    queue: queue_name
  }

  output {
  File seurat_doublet_rds = glob("*.doublet_calls*.rds")[0] 
  }
}



