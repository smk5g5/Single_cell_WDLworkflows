# How to validate the WDL file. 
***This requires that you have the cromwell docker image loaded in an interactive session on compute1.***

***You need to have womtool-53.1.jar in your current working directory if you need to validate/generate json file for inputs. You can download the womtool jar file from here***
https://github.com/broadinstitute/cromwell/releases/download/

I am currently using the following image for cromwell which also happens to be the gms image being used in analysis-workflows pipelines.
```
bsub -Is -q siteman-interactive -G compute-allegra.petti -g /khan.saad/R_seurat -M 128000000 -n 1 -R 'rusage[mem=128000]' -a 'docker(registry.gsc.wustl.edu/apipe-builder/genome_perl_environment:compute1-37)' /bin/bash
```

```
/usr/bin/java -jar womtool-53.1.jar validate ./tasks/single_sample_seurat.wdl
```

# How to write a json file that can then be modified to use with the wdl workflow. 
(***You don't have to write a new json you can modify the existing ones which are there in the example inputs***)

It could be for a single task like the below wdl just does seurat filtering for a single sample

```/usr/bin/java -jar womtool-53.1.jar inputs ./tasks/single_sample_seurat.wdl > test_single_sample_seurat.json```

Or it could be a sub-workflow ![Alt text](./workflow_images/scatter_gather_singleR.png?raw=true "scatter-gather SingleR") e.g. `./subworkflows/scatter_gather_singleR.wdl` which takes a multi-sample/single sample seurat object and runs it against multiple singleR references in a scatter gather fashion and merges their results into the seurat object along with their prediction scores which are stored in a new assay (named based on the input reference name) seurat assay object (singleR references are passed in a tsv as shown in example_inputs file `./example_inputs/SingleR_singleref_scatter_human.tsv` input json looks something like this `./example_inputs/scatter_gather_singleR.json`)

```/usr/bin/java -jar womtool-53.1.jar inputs ./subworkflows/scatter_gather_singleR.wdl > scatter_gather_singleR.json```

or it could be an end to end (multi sample/single sample pipeline,shown here is the multi-sample end to end pipeline) 

![Alt text](./workflow_images/end_to_end_multisample.png?raw=true "End to End multisample workflow")

This multi-sample end-to-end pipeline takes multiple samples,merges them in seurat object ---> Runs doublet calling on each of the sample (10x input) ---> Merges doublet calls for each of the sample in the multi-sample seurat object--->Makes SingleR predictions based on the input singleR references--->Removes the doublet based on majority predictions(ie if majority of the doublet calling methods identify the cell as doublet)--->Renormalizes,reclusters and reruns SingleR to give a final seurat object.

# Running a WDL workflow 
***Here is an example of end to end multi-sample WDL workflow***

***This workflow is run on compute1 as shown here. Change job group, compute-group etc. as necessary.***

```
bsub -oo WDL_end_to_end_multisample_seurat_CT2A.%J.out -G compute-allegra.petti -g /allegrapetti-gms/khan.saad -q siteman -M 8G -R 'select[mem>8G] rusage[mem=8G]' -a 'docker(registry.gsc.wustl.edu/apipe-builder/genome_perl_environment:compute1-37)' /usr/bin/java -Dconfig.file=cromwell_compute1_final.config -jar /opt/cromwell.jar run -t wdl ./pipelines/end_to_end_multisample.wdl -i ./end_to_end_seurat_multisample_CT2A.json
```

# Cromwell-config file

You can modify or use my own cromwell-config that I have here `cromwell_compute1_final.config` 
You will need to change the cromwell logs directory and the root directory as well as the job group(mine is `-g /allegrapetti-gms` you need to change it based on what you see in bjgroup on compute1), compute-group (if you are not using `compute-allegra.petti`). FYI this config file does not have call-caching enabled. Call caching provides the user to be able to restart the WDL workflow from where it failed in case of failure. A caveat for using it is that you need to have a different root directory for every wdl you run since it creates a database lock file which won't be overwritten by a different wdl run and workflow may fail because of that.


Since we are using the gms docker image @chrismiller wrote a script which generates a cromwell-config with call caching enabled that you can use `create_cromwell_config.sh`.

Ideally you would want to run the config generating script from inside the interactive job.
```
bsub -Is -q siteman-interactive -G compute-allegra.petti -g /khan.saad/R_seurat -M 128000000 -n 1 -R 'rusage[mem=128000]' -a 'docker(registry.gsc.wustl.edu/apipe-builder/genome_perl_environment:compute1-37)' /bin/bash
```
Otherwise it would generate paths with `/rdcw` instead of `/storage1` that you would need to replace using sed or from inside vim etc.

Trying to commit code again after dockstore shows up in the settings of the repo
