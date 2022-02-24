version 1.0

import "../tasks/scsorter_run.wdl" as scsorter

workflow scatter_scsorter{

  input {
    String docker_image
    String queue_name
    Int mem_gb
    File scsorter_script
    String Sample_name
    String rds_file_path
  File inputSamplesFile
  Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)
  }

  scatter (sample in inputSamples) {
    call scsorter.run_scsorter as run_scsorter {
      input:
          docker_image=docker_image,
          queue_name=queue_name,
          mem_gb=mem_gb,
          Sample_name=Sample_name,
          rds_file_path=rds_file_path,
          scsorter_script=scsorter_script,
          marker_list_rds=sample[1],
          marker_list_name=sample[0],
  }

}

  call merge_scsorter_results {
      input: 
        docker_image=docker_image, 
        queue_name=queue_name,
        mem_gb=mem_gb,
        Sample_name=Sample_name,
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
    String Sample_name
    String seurat_rds_path
    String type
    Array[File] scsorter_runs
  }
   command <<<
    Rscript ~{merge_scsorter_script} ~{seurat_rds_path} ~{Sample_name} ~{type} ~{sep=" " scsorter_runs}
    >>>
  output {
    File Seurat_merged_scsorter = glob("*.seurat_scsorter_mergedpreds.*.rds")[0]
  }
}

