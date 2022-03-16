version 1.0

task subset_recluster_renormalize{

  input {
    String docker_image
    String queue_name
    Int mem_gb
    File recluster_renormalize_script
    File seurat_rds
    String subset_column_name
    String ident_name
    String inverse
    String output_suffix
  }

   command <<<
    Rscript ~{recluster_renormalize_script} ~{seurat_rds} ~{subset_column_name} ~{ident_name} ~{inverse} ~{output_suffix}
    >>>

  runtime {
    docker : docker_image
    memory: mem_gb + " GB"
    queue: queue_name
  }

  output {
  File seurat_sub_renorm_reclust_rds = glob("*${output_suffix}.rds")[0] 
  }
}

workflow recluster_renormalize{

  call subset_recluster_renormalize
}