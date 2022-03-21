version 1.0

import "../tasks/add_doubletinfo.wdl" as add_doubletinfo
import "../tasks/Doublet_calling.wdl" as Doublet_calling
import "../tasks/single_sample_seurat.wdl" as single_sample_filtering
import "../tasks/clustering_n_pca_simple.wdl" as single_sample_clustering

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



