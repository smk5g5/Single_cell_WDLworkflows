version 1.0

task make_gene_featureplots{

  input {
    String docker_image
    String queue_name
    Int mem_gb
    File gene_feature_plots_script
    File seurat_rds
    File gene_list_file
    String gene_column
    String type_column
  }

   command <<<
    Rscript ~{gene_feature_plots_script} ~{seurat_rds} ~{gene_list_file} ~{gene_column} ~{type_column}
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

workflow Featureplot_genes{

  call make_gene_featureplots
}