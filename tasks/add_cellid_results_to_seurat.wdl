version 1.0

task add_cellid_results_to_seurat{

  input {
    String docker_image
    String queue_name
    Int mem_gb
    File add_cellid_script
    File seurat_rds
    File cellid_predictions
    File cellid_pval_matrix
    String marker_file_reference_name
  }

   command <<<
    Rscript ~{add_cellid_script} ~{seurat_rds} ~{cellid_predictions} ~{cellid_pval_matrix} ~{marker_file_reference_name}
    >>>

  runtime {
    docker : docker_image
    memory: mem_gb + " GB"
    queue: queue_name
  }

  output {
 File output_rds = glob("*.Cellid_results.*.rds")[0]
  }
}

workflow Add_cellid_to_srt{

  call add_cellid_results_to_seurat
}