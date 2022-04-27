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
    String organism
  }

   command <<<
    Rscript ~{recluster_renormalize_script} ~{seurat_rds} ~{subset_column_name} ~{ident_name} ~{inverse} ~{output_suffix} ~{organism}
    >>>

  runtime {
    docker : docker_image
    memory: mem_gb + " GB"
    queue: queue_name
  }

  output {
  File seurat_sub_renorm_reclust_rds = glob("*${output_suffix}*.rds")[0] 
  }
}

workflow recluster_renormalize{

  call subset_recluster_renormalize
}

version 1.0

workflow linear_chain_single_sample_seurat_doubletsadd{

    input {
    String docker_image
    String queue_name
    Int mem_gb
    String organism    
    String Sample_name
    String cellranger_outs_directory
    }
    call Doublet_calling.run_doublet_collection as run_doublet {
    input:
    docker_image=docker_image,
    queue_name=queue_name,
    mem_gb=mem_gb,
    cellranger_outs_directory=cellranger_outs_directory,
    Sample_name=Sample_name
    }
    call single_sample_filtering.run_seurat_singlesample as run_single_srt {
    input:
    docker_image=docker_image,
    queue_name=queue_name,
    mem_gb=mem_gb,
    cellranger_outs_directory=cellranger_outs_directory,
    Sample_name=Sample_name
    }
    call single_sample_clustering.run_clustering_n_pca_simple as run_srt_clust_simp {
    input:
    docker_image=docker_image,
    queue_name=queue_name,
    mem_gb=mem_gb,
    rds_file_path=run_single_srt.intermed_rds,
    Sample_name=Sample_name
    }
    call add_doubletinfo.add_doublets_metadata as add_doub_to_srt {
    input:
    docker_image=docker_image,
    queue_name=queue_name,
    mem_gb=mem_gb,
    input_rds_file=run_srt_clust_simp.clustering_rds,
    doublet_file=run_doublet.doublet_results
    }
    output {
        File doubmeta_rds=add_doub_to_srt.seurat_doublet_rds
    }
}