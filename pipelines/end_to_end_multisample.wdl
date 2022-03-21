version 1.0

import "../tasks/multisample_seurat.wdl" as multisample_seurat
import "../subworkflows/scatter_gather_Doubletcall_multisample.wdl" as scatter_gather_doublet
import "../subworkflows/scatter_gather_singleR.wdl" as scatter_gather_singleR
import "../subworkflows/recluster_renorm_rerun_singleR.wdl" as recluster_renorm_rerun_singleR

workflow end_to_end_seurat_single_sample{

    input {
    String docker_image
    String queue_name
    Int mem_gb
    File inputSamplesFile
    File singleR_refsFile
    String subset_column_name
    String ident_name
    String inverse
    String output_suffix
    String organism    
    File recluster_renormalize_script
    Array[Array[File]] inputSamples = read_tsv(inputSamplesFile)
    }

    call multisample_seurat.run_seurat_multisample as run_seurat_multisample {
        input:
            docker_image=docker_image,
            queue_name=queue_name,
            mem_gb=mem_gb,
            organism=organism,
            project_name=project_name,
            gene_lists=gene_lists,
            seurat_multisample_rscript=seurat_multisample_rscript,
            multisample_seurat_10x_inp=inputSamplesFile,
    }
    call scatter_gather_doublet.scatter_doublet as scatter_gather_doublet {
        input:
            docker_image=docker_image,
            queue_name=queue_name,
            mem_gb=mem_gb,
            inputSamplesFile=inputSamplesFile,
            merge_doublet_calls_in_seurat_script=merge_doublet_calls_in_seurat_script,
            multisample_seurat_rds=run_seurat_multisample.final_rds,
    }
    call scatter_gather_singleR.scatter_gather_singleR as scatter_gather_singleR {
        input:
            docker_image=docker_image,
            queue_name=queue_name,
            mem_gb=mem_gb,
            inputSamplesFile=inputSamplesFile,
            singleR_singleref_rscript=singleR_singleref_rscript,
            merge_doublet_calls_in_seurat_script=merge_doublet_calls_in_seurat_script,
            multisample_seurat_rds=scatter_gather_doublet.seurat_doublet_rds,
    }
    call recluster_renorm_rerun_singleR.LinearChain_recluster_rerun_singleR as LinearChain_recluster_rerun_singleR {
        input:
        docker_image=docker_image,
        queue_name=queue_name,
        mem_gb=mem_gb,
        seurat_rds=scat_gath_singleR.seurat_singleR_rds,
        subset_column_name=subset_column_name,
        ident_name=ident_name,
        inverse=inverse,
        organism=organism,
        output_suffix=output_suffix,
        inputSamplesFile=singleR_refsFile
        }

}