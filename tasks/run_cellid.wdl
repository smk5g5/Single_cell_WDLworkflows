version 1.0

task run_cellid{

  input {
    String docker_image
    String queue_name
    Int mem_gb
    File cellid_run_rscript
    File seurat_rds
    File marker_gene_list_rds
    String marker_file_reference_name
  }

   command <<<
    Rscript ~{cellid_run_rscript} ~{seurat_rds} ~{marker_gene_list_rds} ~{marker_file_reference_name}
    >>>

  runtime {
    docker : docker_image
    memory: mem_gb + " GB"
    queue: queue_name
  }

  output {
 File cellid_predictions = glob("CellID_final_predictions.*.rds")[0]
 File cellid_pval_matrix = glob("CellID_pval_matrix.*.rds")[0]
  }
}

workflow CELLID_run{

  call run_cellid
}