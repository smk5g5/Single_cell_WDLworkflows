import scvi
import scanpy as sc
import pandas as pd
import seaborn as sns
import numpy as np
import os 
os.environ['MPLCONFIGDIR'] = os.getcwd() + "/configs/"
#this will make the config directory to be the PWD instead of $HOME
import matplotlib.pyplot as plt
import numpy as np
from scipy.io import mmread
import scipy.sparse as sp
import argparse

# cell_cycle_genes = [x.strip() for x in open('./data/regev_lab_cell_cycle_genes.txt')]


#This script will be updated to work with WDL
#for smk5g5/scvi-tools:latest call python3 not python

if __name__ == "__main__":
	parser=argparse.ArgumentParser(description="""
	This script will run cellassign on counts data geneXcell
    matrix with genes as rows and cells as columns""")

	parser.add_argument("Counts", help="Count input file")
	parser.add_argument("sample_name", help="sample name")
	parser.add_argument("organism", help="organism name")

	args = parser.parse_args()

	COUNTS=args.Counts
	Sample_name=args.sample_name
	Organism=args.organism
	output_h5ad='anndata_'+Sample_name+'.h5ad'
	f_exprMat = str(COUNTS)
	print("f_exprMat ",f_exprMat)
	adata = sc.read_text(f_exprMat, delimiter='\t', first_column_names=True )
	# adata = adata.transpose()
	#print top 20 genes with highest expression
	sc.pl.highest_expr_genes(adata, n_top=20,save=True)
	#Not filtering with already filtered seurat object
	# sc.pp.filter_cells(adata, min_genes=200)
	# sc.pp.filter_genes(adata, min_cells=3)

	#Not removing mitochondrial genes as well (since its a filtered seurat object)
	# adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
	# sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
	adata.var['mt'] = adata.var_names.str.startswith(('MT-', 'mt-'))  # annotate the group of mitochondrial genes as 'mt'
	sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
	fig_out = '_'+Sample_name+'.pdf'  #append sample name at the end of plot name
	sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],jitter=0.4, multi_panel=True,save=fig_out)
	fig_out = '_total_counts_vs_pct_counts_mt_'+Sample_name+'.pdf'
	sc.pl.scatter(adata, x='total_counts', y='pct_counts_mt',save=fig_out)
	fig_out = '_total_counts_vs_n_genes_by_counts_'+Sample_name+'.pdf'
	sc.pl.scatter(adata, x='total_counts', y='n_genes_by_counts',save=fig_out)
	#log normalize the data
	sc.pp.normalize_total(adata, target_sum=1e4)
	sc.pp.log1p(adata)
	#identify highly variable genes
	sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
	#plot highly variable genes
	sc.pl.highly_variable_genes(adata)
	adata.raw = adata
	adata.var['total_expression'] = adata.X.sum(0)

	if(str(Organism))=='human':
		cell_cycle_df = pd.read_csv("/storage1/fs1/allegra.petti/Active/10xGenomics/key.gene.lists/CellCycleTirosh.txt",sep="\t",header=0)
		s_genes = list(cell_cycle_df['Gene.Symbol'][cell_cycle_df['List']=='G1/S'])
		g2m_genes = list(cell_cycle_df['Gene.Symbol'][cell_cycle_df['List']=='G2/M'])
		cell_cycle_genes = list(cell_cycle_df['Gene.Symbol'])
		cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
		# cell_cycle_genes = [x.strip() for x in open('/data/regev_lab_cell_cycle_genes.txt')]
	else:
		cell_cycle_df = pd.read_csv("/storage1/fs1/allegra.petti/Active/10xGenomics/key.gene.lists/CellCycleTirosh_mouse.txt",sep="\t",header=None)
		cell_cycle_df =['List', 'Gene.Symbol']
		s_genes = list(cell_cycle_df['Gene.Symbol'][cell_cycle_df['List']=='G1/S'])
		g2m_genes = list(cell_cycle_df['Gene.Symbol'][cell_cycle_df['List']=='G2/M'])

	sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
	# Regress out effects of cellcycle
	sc.pp.regress_out(adata, ['S_score', 'G2M_score'],n_jobs=10)
	sc.pp.scale(adata)
	# Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed. Scale the data to unit variance.
	# sc.pp.regress_out(adata, ['total_counts', 'pct_counts_mt']) #might not need to do this
	# sc.pp.scale(adata, max_value=10)
	sc.tl.pca(adata, svd_solver='arpack')
	fig_out = "_"+Sample_name+'.pdf'
	sc.pl.pca_variance_ratio(adata, log=True,save=fig_out)
	sc.pp.neighbors(adata, n_neighbors=10, n_pcs=30)
	top10_genes = adata.var['total_expression'].nlargest(10)
	top10_genes=list(top10_genes.index)
	sc.tl.umap(adata)
	fig_out = '_top10_genes_'+Sample_name+'.pdf'

	sc.pl.umap(adata, color=top10_genes,save=fig_out)

	sc.tl.leiden(adata)
	fig_out = 'leiden_'+Sample_name+'.pdf'

	sc.pl.umap(adata, color=['leiden'],save=fig_out)

	adata.write(output_h5ad)