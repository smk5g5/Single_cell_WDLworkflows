version 1.0

import "../tasks/SingleR_singleref.wdl" as singleR

workflow scatter_gather_singleR_singleref{

  input {
    String docker_image
    String queue_name
    Int mem_gb
    File singleR_singleref_rscript
    String singleR_ref_rds
    String reference_name
    String label_column_name
    File singleR_tsv_file
  File inputSamplesFile
  Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)
  }

  scatter (sample in inputSamples) {
    call singleR.run_singleR_singleref as run_singleR {
      input:
          docker_image=docker_image,
          queue_name=queue_name,
          mem_gb=mem_gb,
          seurat_rds=sample[1],
          singleR_singleref_rscript=singleR_singleref_rscript,
          reference_name=reference_name,
          label_column_name=label_column_name,
          singleR_ref_rds = singleR_ref_rds
  }
    call add_singleR_results_to_seurat {
        input:
          docker_image=docker_image,
          queue_name=queue_name,
          mem_gb=mem_gb,
          merge_SingleR_in_seurat_script=merge_SingleR_in_seurat_script,
          singleR_tsv_file=singleR_tsv_file,
          seurat_rds=sample[1],
          singleR_pred_rds=run_singleR.singleR_pred_rds,
    }
  output {
  File seurat_singleR_rds = add_singleR_results_to_seurat.seurat_singleR_rds
  }  
}
}

task add_singleR_results_to_seurat{

  input {
    String docker_image
    String queue_name
    Int mem_gb
    File merge_SingleR_in_seurat_script
    File seurat_rds
    File singleR_tsv_file
    Array[String] singleR_pred_rds
  }

   command <<<
    Rscript ~{merge_SingleR_in_seurat_script} ~{seurat_rds} ~{singleR_tsv_file} ~{sep=" " singleR_pred_rds}
    >>>

  runtime {
    docker : docker_image
    memory: mem_gb + " GB"
    queue: queue_name
  }

  output {
  File seurat_singleR_rds = glob("*.single_R_preds.*.rds")[0] 
  }
}
