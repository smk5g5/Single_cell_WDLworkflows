# end_to_end_multisample.wdl 

![Alt text](./workflow_images/end_to_end_multisample.png?raw=true "End to End multisample workflow")

This multi-sample end-to-end pipeline takes multiple samples,merges them in seurat object ---> Runs doublet calling on each of the sample (10x input) ---> Merges doublet calls for each of the sample in the multi-sample seurat object--->Makes SingleR predictions based on the input singleR references--->Removes the doublet based on majority predictions(ie if majority of the doublet calling methods identify the cell as doublet)--->Renormalizes,reclusters and reruns SingleR to give a final seurat object.

# end_to_end_seurat_singlesample.wdl 

![Alt text](./workflow_images/end_to_end_seurat_singlesample.png?raw=true "End to End single sample seurat workflow")

Unlike multi-sample seurat workflows this makes SingleR calls on individual samples --> removes doublets from each sample --> recalls SingleR --> Merges the single sample into one multi-sample seurat object and renormalizes and reclusters giving a final seurat object.
