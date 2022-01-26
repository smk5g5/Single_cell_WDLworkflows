#How to validate the WDL file. This requires that you have the cromwell docker image loaded in an interactive session on compute1.

```
/usr/bin/java -jar /scratch1/fs1/allegra.petti/khan.saad/WDL_workflow/womtool-53.1.jar validate ./tasks/single_sample_seurat.wdl
```

#How to write a json file that can then be modified to use with the wdl workflow. 

```/usr/bin/java -jar /scratch1/fs1/allegra.petti/khan.saad/WDL_workflow/womtool-53.1.jar inputs ./tasks/single_sample_seurat.wdl > test_single_sample_seurat.json```


/usr/bin/java -jar /scratch1/fs1/allegra.petti/khan.saad/WDL_workflow/womtool-53.1.jar inputs ./tasks/clustering_n_pca_simple.wdl > test_run_clustering_n_pca_simple.json


#Running a WDL workflow interactively
```
bsub -Is -G compute-bolton -g /bwileytest -q general-interactive -M 8G -R 'select[mem>8G] rusage[mem=8G]' -a 'docker(broadinstitute/cromwell:dev)' /bin/bash

/opt/java/openjdk/bin/java \
    -Dconfig.file=/storage1/fs1/bolton/Active/data/presets/cromwell.config \
    -jar /app/cromwell.jar run \
    -o /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/example_data/Brian/options/vardict_options.json \
    -t wdl \
    -i /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/example_data/tools/vardict_tool.json \
    /storage1/fs1/bolton/Active/projects/TwinStrand/PRJ00087.2021.05.24.deliverables/wdl/analysis-wdls/definitions/subworkflows/vardict.wdl
```

#Running a WDL workflow non-interactively

```
bsub -oo test_WDL_workflows_simpleclus.%J.out -G compute-allegra.petti -g /allegrapetti-gms/khan.saad -q general -M 8G -R 'select[mem>8G] rusage[mem=8G]' -a 'docker(registry.gsc.wustl.edu/apipe-builder/genome_perl_environment:compute1-37)' /usr/bin/java -Dconfig.file=/scratch1/fs1/allegra.petti/khan.saad/WDL_workflow/Single_cell_WDLworkflows/cromwell_compute1.config -jar /opt/cromwell.jar run -t wdl ./tasks/clustering_n_pca_simple.wdl -i ./test_run_clustering_n_pca_simple.json
```
/storage1/fs1/allegra.petti/Active/Users/khan.saad/WDL_pipelines/Seurat_single_sample/04b5857b-e769-4b63-95c9-734067908c6c/call-run_seurat_singlesample/execution/B148.MergedFilteredSeuratObject.20220112.rds


```
/scratch1/fs1/allegra.petti/khan.saad/WDL_workflow/Single_cell_WDLworkflows/test_run_clustering_n_pca_simple.json
```