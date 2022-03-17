version 1.0

import "../tasks/subset_renormalize_recluster.wdl" as recluster_renormalize

workflow LinearChain_recluster_rerun_singleR{

  input {
  String docker_image
    String queue_name
    Int mem_gb
    File singleR_singleref_rscript
    File merge_SingleR_in_seurat_script
    File seurat_rds
  File inputSamplesFile
  Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)
  }

  call recluster_renormalize.subset_recluster_renormalize { input: seurat_rds=seurat_rds }

  scatter (sample in inputSamples) {
    call singleR.run_singleR_singleref as run_singleR {
      input:
          docker_image=docker_image,
          queue_name=queue_name,
          mem_gb=mem_gb,
          singleR_singleref_rscript=singleR_singleref_rscript,
          seurat_rds=recluster_renormalize.subset_recluster_renormalize.seurat_sub_renorm_reclust_rds,
          reference_name=sample[0],
          singleR_ref_rds=sample[2],
        label_column_name=sample[1]
  }}

 call add_singleR_results_to_seurat {
      input: 
        docker_image=docker_image, 
        queue_name=queue_name,
        mem_gb=mem_gb,
        merge_SingleR_in_seurat_script=merge_SingleR_in_seurat_script,
        inputSamplesFile=inputSamplesFile,
        seurat_rds=recluster_renormalize.subset_recluster_renormalize.seurat_sub_renorm_reclust_rds,
        singleR_pred_rds=run_singleR.singleR_pred_rds
    }
}


task add_singleR_results_to_seurat{

  input {
    String docker_image
    String queue_name
    Int mem_gb
    File merge_SingleR_in_seurat_script
    File inputSamplesFile
    File seurat_rds
    Array[String] singleR_pred_rds
  }

   command <<<
    Rscript ~{merge_SingleR_in_seurat_script} ~{seurat_rds} ~{inputSamplesFile} ~{sep=" " singleR_pred_rds}
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