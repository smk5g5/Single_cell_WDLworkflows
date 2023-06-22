version 1.0

task run_singleR_singleref {

  input {
    String docker_image
    String queue_name
    Int mem_gb
    File singleR_singleref_rscript
    File seurat_rds
    File singleR_ref_rds
    String reference_name
    String label_column_name
    String Sample_name
  }

   command <<<
    Rscript ~{singleR_singleref_rscript} ~{seurat_rds} ~{singleR_ref_rds} ~{reference_name} ~{label_column_name} ~{Sample_name}
    >>>

  runtime {
    docker : docker_image
    memory: mem_gb + " GB"
    queue: queue_name
  }

  output {
  File singleR_pred_rds = glob("*_singleR_preds.rds")[0]
  }
}

workflow singleR_singleref {

  call run_singleR_singleref
}
