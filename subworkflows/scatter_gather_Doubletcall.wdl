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

        seurat_rds_path=rds_file_path,
        scsorter_runs=run_scsorter.scsorter_preds_rds
    }

}

task merge_scsorter_results {
  input {
    String docker_image
    String queue_name
    Int mem_gb
    File merge_scsorter_script
    File inputSamplesFile
    String project_name
    String seurat_rds_path
    String split_by
    Array[String] scsorter_runs
  }
   command <<<
    Rscript ~{merge_scsorter_script} ~{seurat_rds_path} ~{project_name} ~{inputSamplesFile} ~{split_by} ~{sep=" " scsorter_runs}
    >>>

  runtime {
      docker : docker_image
      memory: mem_gb + " GB"
      queue: queue_name
  }

  output {
    File Seurat_merged_scsorter = glob("*.seurat_scsorter_mergedpreds.*.rds")[0]
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
  File seurat_doublet_rds = glob("*.doublet_calls.rds")[0] 
  }
}



