version 1.0

task run_doublet_collection{

  input {
    String docker_image
    String queue_name
    Int mem_gb
    File doublet_calling_script
    String cellranger_outs_directory
    String Sample_name
  }

   command <<<
    Rscript ~{doublet_calling_script} ~{cellranger_outs_directory} ~{Sample_name}
    >>>

  runtime {
    docker : docker_image
    memory: mem_gb + " GB"
    queue: queue_name
  }

  output {
  File doublet_results = glob("*_Doublet_collection_results.rds")[0] 
  }
}

workflow Doublet_colletion{

  call run_doublet_collection
}