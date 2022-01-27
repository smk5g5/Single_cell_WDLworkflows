version 1.0

task run_singleR_scsorter {

  input {
    String docker_image
    String queue_name
    Int mem_gb
    File singleR_scsorter_rscript
    String rds_file_path
    String Sample_name
    Float clustering_res
    Int nPC
    String singleR_yaml
  }

   command <<<
    Rscript ~{singleR_scsorter_rscript} ~{rds_file_path} ~{Sample_name} ~{clustering_res} ~{nPC} ~{singleR_yaml}
    >>>

  runtime {
    docker : docker_image
    memory: mem_gb + " GB"
    queue: queue_name
  }

  output {
  Array[File] output_rds = glob("*.rds")
  Array[File] output_jpgs = glob("*.jpg")
  File intermed_rds = glob("*.MergedFilteredSeuratObject.*.rds")[0] 
  }
}

workflow SingleR_scsorter {

  call run_singleR_scsorter
}