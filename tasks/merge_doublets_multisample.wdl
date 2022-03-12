version 1.0

task run_seurat_multisample{

  input {
    String docker_image
    String queue_name
    Int mem_gb
    File multisample_seurat_10x_inp
    String organism
    String project_name
    File gene_lists
    File seurat_multisample_rscript
  }
   command <<<
    Rscript ~{seurat_multisample_rscript} ~{multisample_seurat_10x_inp} ~{organism} ~{project_name} ~{gene_lists}
    >>>

  runtime {
    docker : docker_image
    memory: mem_gb + " GB"
    queue: queue_name
  }

  output {
  Array[File] output_pdfs = glob("*.pdf")
  Array[File] output_pngs = glob("*.png")
  Array[File] output_rds = glob("*.rds")
  Array[File] output_jpgs = glob("*.jpg")
  File final_rds = glob("*.SCT.PCA.UMAP.TSNE.CLUST.*.rds")[0] 
  }
}

workflow Seurat_multisample {

  call run_seurat_multisample
}

