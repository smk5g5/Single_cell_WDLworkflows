import scvi
import scanpy as sc
import pandas as pd
import seaborn as sns
import numpy as np
import os os.environ['MPLCONFIGDIR'] = os.getcwd() + "/configs/"
#this will make the config directory to be the PWD instead of $HOME
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import mmread
import scipy.sparse as sp

import os
import argparse

#This script will be updated
#for smk5g5/scvi-tools:latest call python3 not python

if __name__ == "__main__":
	parser=argparse.ArgumentParser(description="""
	This script will run cellassign on counts data geneXcell
    matrix with genes as rows and cells as columns""")

	parser.add_argument("Counts", help="Count input file")
	parser.add_argument("sample_name", help="sample name")

	args = parser.parse_args()

	COUNTS=args.Counts
	Sample_name=args.sample_name
	output_h5ad='anndata_'+Sample_name+'.h5ad'
	f_exprMat = str(COUNTS)
	print("f_exprMat ",f_exprMat)
	adata = sc.read_text(f_exprMat, delimiter='\t', first_column_names=True )
	adata = adata.transpose()
	#print top 20 genes with highest expression
	sc.pl.highest_expr_genes(adata, n_top=20, )
	#Not filtering with already filtered seurat object
	# sc.pp.filter_cells(adata, min_genes=200)
	# sc.pp.filter_genes(adata, min_cells=3)

	#Not removing mitochondrial genes as well (since its a filtered seurat object)
	# adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
	# sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
	adata.var['mt'] = adata.var_names.str.startswith(('MT-', 'mt-'))  # annotate the group of mitochondrial genes as 'mt'
	sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
	fig_out = 'violin_plot_'+Sample_name+'.pdf'
	sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],jitter=0.4, multi_panel=True,save=fig_out)
	fig_out = 'scatter_plot_'+Sample_name+'.pdf'
	sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt',save=fig_out)
	fig_out = 'scatter_plot_'+Sample_name+'.pdf'
	sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',save=fig_out)
	#log normalize the data
	sc.pp.normalize_total(adata, target_sum=1e4)
	sc.pp.log1p(adata)
	#identify highly variable genes
	sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
	#plot highly variable genes
	sc.pl.highly_variable_genes(adata)
	adata.raw = adata

	# Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed. Scale the data to unit variance.
	# sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt']) #might not need to do this
	sc.pp.scale(adata, max_value=10)
	sc.tl.pca(adata, svd_solver='arpack')
	fig_out = 'pca_variance_ratio_'+Sample_name+'.pdf'
	sc.pl.pca_variance_ratio(adata, log=True,save=fig_out)
	sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)

	top10_genes = adata.var['total_expression'].nlargest(10)
	sc.tl.umap(adata)
	fig_out = 'top10_genes_'+Sample_name+'.pdf'

	sc.pl.umap(adata, color=top10_genes,save=fig_out)

	sc.tl.leiden(adata)
	fig_out = 'leiden_clustering_'+Sample_name+'.pdf'

	sc.pl.umap(adata, color=['leiden', top10_genes],save=fig_out)

	adata.write(output_h5ad)