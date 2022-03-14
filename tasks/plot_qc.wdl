version 1.0

task plot_qc{

  input {
    String docker_image
    String queue_name
    Int mem_gb
    File plot_qc_script
    File seurat_rds
    String organism
    String project_name
    String pre_post
  }

   command <<<
    Rscript ~{plot_qc_script} ~{seurat_rds} ~{organism} ~{project_name} ~{pre_post}
    >>>

  runtime {
    docker : docker_image
    memory: mem_gb + " GB"
    queue: queue_name
  }

  output {
  Array[File] output_pdfs = glob("*.pdf")
  }
}

workflow QC_plots{

  call plot_qc
}