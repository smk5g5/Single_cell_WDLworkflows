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
  }

   command <<<
    Rscript ~{htmap_rscript} ~{Seurat_rds} ~{DEG_file} ~{output_prefix}
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

workflow DEGs_downstream {

  call make_heatmap_DEGs
}