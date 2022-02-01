version 1.0

task make_heatmap_DEGs{

  input {
    String docker_image
    String queue_name
    Int mem_gb
    File htmap_rscript
    String Seurat_rds
    String output_prefix
    String DEG_file
    String clustering
    String aggregate_exp_rds
  }

   command <<<
    Rscript ~{htmap_rscript} ~{Seurat_rds} ~{DEG_file} ~{output_prefix} ~{clustering} ~{aggregate_exp_rds}
    >>>

  runtime {
    docker : docker_image
    memory: mem_gb + " GB"
    queue: queue_name
  }

  output {
  Array[File] output_jpgs = glob("*.jpg")
  Array[File] output_rds = glob("*.rds")
  Array[File] output_pngs = glob("*.png")
  }
}


task aggregate_expression_seurat4_0_3{

  input {
    String docker_image
    String queue_name
    Int mem_gb
    File aggregate_script
    String Seurat_rds
    String output_prefix
    String DEG_file
    String clustering
  }

   command <<<
    Rscript ~{aggregate_script} ~{Seurat_rds} ~{clustering} ~{DEG_file} 
    >>>

  runtime {
    docker : docker_image
    memory: mem_gb + " GB"
    queue: queue_name
  }

  output {
  Array[File] aggregate_rds = glob("aggregate_expression_obj.rds")[0]
  }
}

workflow DEGs_downstream {
  call aggregate_expression_seurat4_0_3 
  call make_heatmap_DEGs{
    input: 
    Seurat_rds = aggregate_expression_seurat4_0_3.Seurat_rds,
    clustering = aggregate_expression_seurat4_0_3.clustering,
    DEG_file = aggregate_expression_seurat4_0_3.DEG_file,
    aggregate_exp_rds = aggregate_expression_seurat4_0_3.aggregate_exp_rds
  }
}
