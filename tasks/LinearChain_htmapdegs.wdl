workflow DEGs_downstream {
  String Seurat_rds
  String output_prefix
  String DEG_file
  String clustering
  call aggregate_expression_seurat4_0_3 { 
    input: 
    Seurat_rds=Seurat_rds,
    DEG_file=DEG_file,
    clustering=clustering}
  call make_heatmap_DEGs {
   input: 
   Seurat_rds=Seurat_rds,
   output_prefix=output_prefix,
   DEG_file=DEG_file,
   clustering=clustering,
   aggregate_exp_rds=aggregate_expression_seurat4_0_3.aggregate_rds
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
