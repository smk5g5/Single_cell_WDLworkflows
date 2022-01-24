version 1.0

task run_seurat_singlesample{

  input {
    String docker_image
    String queue_name
    Int mem_gb
    File seurat_singlesample_rscript
    String cellranger_outs_directory
    String Sample_name
  }

   command <<<
    Rscript ~{seurat_singlesample_rscript} ~{cellranger_outs_directory} ~{Sample_name}
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
  Array[File] intermed_rds = glob("*.MergedFilteredSeuratObject.*.rds") 
  }

}

workflow Seurat_single_sample {

  call run_seurat_singlesample
}

