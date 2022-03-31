# How to validate the WDL file. This requires that you have the cromwell docker image loaded in an interactive session on compute1.

I am currently using the following image for cromwell which also happens to be the latest gms image
```
bsub -Is -q siteman-interactive -G compute-allegra.petti -g /khan.saad/R_seurat -M 128000000 -n 1 -R 'rusage[mem=128000]' -a 'docker(registry.gsc.wustl.edu/apipe-builder/genome_perl_environment:compute1-37)' /bin/bash
```

```
/usr/bin/java -jar /scratch1/fs1/allegra.petti/khan.saad/WDL_workflow/womtool-53.1.jar validate ./tasks/single_sample_seurat.wdl
```

# How to write a json file that can then be modified to use with the wdl workflow. 

```/usr/bin/java -jar /scratch1/fs1/allegra.petti/khan.saad/WDL_workflow/womtool-53.1.jar inputs ./tasks/single_sample_seurat.wdl > test_single_sample_seurat.json```


/usr/bin/java -jar /scratch1/fs1/allegra.petti/khan.saad/WDL_workflow/womtool-53.1.jar inputs ./tasks/clustering_n_pca_simple.wdl > test_run_clustering_n_pca_simple.json


# Running a WDL workflow 
# Here is an example of end to end multi-sample WDL workflow
![Alt text](./workflow_images/end_to_end_multisample.png?raw=true "End to End multisample workflow")
# This workflow is run on compute1 as shown here.

```
bsub -oo WDL_end_to_end_multisample_seurat_CT2A.%J.out -G compute-allegra.petti -g /allegrapetti-gms/khan.saad -q siteman -M 8G -R 'select[mem>8G] rusage[mem=8G]' -a 'docker(registry.gsc.wustl.edu/apipe-builder/genome_perl_environment:compute1-37)' /usr/bin/java -Dconfig.file=/scratch1/fs1/allegra.petti/khan.saad/WDL_workflow/Single_cell_WDLworkflows/cromwell_compute1_final.config -jar /opt/cromwell.jar run -t wdl ./pipelines/end_to_end_multisample.wdl -i ./end_to_end_seurat_multisample_CT2A.json
```
