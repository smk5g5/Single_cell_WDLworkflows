#How to validate the WDL file. This requires that you have the cromwell docker image loaded in an interactive session on compute1.

```
/usr/bin/java -jar /scratch1/fs1/allegra.petti/khan.saad/WDL_workflow/womtool-53.1.jar validate ./tasks/single_sample_seurat.wdl
```

#How to write a json file that can then be modified to use with the wdl workflow. 

```/usr/bin/java -jar /scratch1/fs1/allegra.petti/khan.saad/WDL_workflow/womtool-53.1.jar inputs ./tasks/single_sample_seurat.wdl > test_single_sample_seurat.json```