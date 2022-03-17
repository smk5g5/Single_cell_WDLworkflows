version 1.0

import "../subworkflows/scatter_gather_singleR.wdl" as scatter_gather_singleR
import "../tasks/subset_renormalize_recluster.wdl" as recluster_renormalize

workflow LinearChain_recluster_rerun_singleR {
    input {
    String docker_image
    String queue_name
    Int mem_gb
    File seurat_rds
}
  call recluster_renormalize.subset_recluster_renormalize { input: seurat_rds=seurat_rds }
  call scatter_gather_singleR.run_singleR { input: seurat_rds=recluster_renormalize.subset_recluster_renormalize.seurat_sub_renorm_reclust_rds }
  call scatter_gather_singleR.add_singleR_results_to_seurat {
   input: 
   seurat_rds=recluster_renormalize.subset_recluster_renormalize.seurat_sub_renorm_reclust_rds,
   singleR_pred_rds=scatter_gather_singleR.run_singleR.singleR_pred_rds
  }
}