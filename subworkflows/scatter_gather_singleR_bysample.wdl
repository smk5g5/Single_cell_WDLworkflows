version 1.0

import "../tasks/SingleR_singleref.wdl" as singleR_run

workflow run_singleR_by_sample{
File filtered_vcfs_list
Array[File] input_vcfs = read_lines(filtered_vcfs_list)
  input {
  String docker_image
    String queue_name
    Int mem_gb
    File singleR_singleref_rscript
    File merge_SingleR_in_seurat_script
    File singleR_ref_rds
    String reference_name
    String label_column_name
    File singleR_tsv_file
    File inputSamplesFile
   Array[File] inputSamples = read_lines(inputSamplesFile)
  }

  scatter (sample in inputSamples) {
    call singleR_run.run_singleR_singleref as run_singleR {
      input:
      docker_image=docker_image,
      queue_name=queue_name,
      mem_gb=mem_gb,
      singleR_singleref_rscript=singleR_singleref_rscript,
      seurat_rds=sample,
      reference_name=reference_name,
      label_column_name=label_column_name,
      singleR_ref_rds=singleR_ref_rds
    }
    call add_singleR_results_to_seurat {
      input: 
        docker_image=docker_image, 
        queue_name=queue_name,
        mem_gb=mem_gb,
        merge_SingleR_in_seurat_script=merge_SingleR_in_seurat_script,
        singleR_tsv_file=singleR_tsv_file,
        seurat_rds=sample,
        singleR_pred_rds=run_singleR.singleR_pred_rds
    }

  }
}

task add_singleR_results_to_seurat{

  input {
    String docker_image
    String queue_name
    Int mem_gb
    File merge_SingleR_in_seurat_script
    File singleR_tsv_file
    File seurat_rds
    File singleR_pred_rds
  }

   command <<<
    Rscript ~{merge_SingleR_in_seurat_script} ~{seurat_rds} ~{singleR_tsv_file} ~{singleR_pred_rds}
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