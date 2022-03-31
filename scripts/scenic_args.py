import os
import numpy as np
import pandas as pd
import scanpy as sc
import loompy as lp
from MulticoreTSNE import MulticoreTSNE as TSNE
import umap
from MulticoreTSNE import MulticoreTSNE as TSNE
import importlib
import seaborn as sns
import matplotlib.pyplot as plt

#This script will be updated to work with WDL

import argparse

if __name__ == "__main__":
	parser=argparse.ArgumentParser(description="""
	This script will run Stream bifurcation analysis
	pass Counts file metadata file output dir and sample name
	""")
	parser.add_argument("Counts", help="Count input file")
	parser.add_argument("meta", help="metadata input file")
	parser.add_argument("outdir", help="output directory")
	parser.add_argument("metaname", help="output_file_meta name")

	args = parser.parse_args()

	COUNTS=args.Counts
	META=args.meta
	OUTDIR=args.outdir
	metaname=args.metaname


	print("transposed Counts input : " + COUNTS)
	print("metadata input : " + META)
	print("Output directory : " + OUTDIR)
	print("metaname output: " + metaname)

	sc.settings.njobs = 20

	wdir = OUTDIR

	os.chdir( wdir )

	f_exprMat = str(wdir+'/'+COUNTS)
	f_metadata = str(wdir+'/'+META)
	f_loom_path_unfilt = str(wdir+'/'+metaname+'.final.unfilt.loom')
	f_loom_path_scenic = str(wdir+'/'+metaname+'.final.reanalyzed.loom') # schwann data reanalyzed with scanpy
	f_anndata_path = str(wdir+'/'+metaname+'.anndata.h5ad')
	f_pyscenic_output = str(wdir+'/'+metaname+'.scenic.loom')
	f_final_loom = str(wdir+'/'+metaname+'.scenic.scope.loom')

	print("f_exprMat ",f_exprMat)
	print("f_metadata ",f_metadata)
	print("f_loom_path_unfilt ",f_loom_path_unfilt)
	print("f_loom_path_scenic ",f_loom_path_scenic)
	print("f_anndata_path ",f_anndata_path)
	print("f_pyscenic_output ",f_pyscenic_output)
	print("f_final_loom ",f_final_loom)

	adata = sc.read_text( f_exprMat, delimiter='\t', first_column_names=True )

	row_attrs = {
	"Gene": np.array(adata.var.index) ,
	}
	col_attrs = {
	"CellID":  np.array(adata.obs.index) ,
	"nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
	"nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
	}

	# this creates the loom file containing the unfiltered data:
	lp.create( f_loom_path_unfilt, adata.X.transpose(), row_attrs, col_attrs ) # creating a loom scrna object here, but for unfiltered data

	### Initial/basic filtering
	# [6] now read it in and process it:

	adata = sc.read_loom( f_loom_path_unfilt )

	nCountsPerGene = np.sum(adata.X, axis=0)
	nCellsPerGene = np.sum(adata.X>0, axis=0)

	print("Number of counts (in the dataset units) per gene:", nCountsPerGene.min(), " - " ,nCountsPerGene.max())
	print("Number of cells in which each gene is detected:", nCellsPerGene.min(), " - " ,nCellsPerGene.max())

	nCells = adata.X.shape[0]

	minCountsPerGene=3*.01*nCells # 3 counts in 1% of cells
	print("minCountsPerGene: ", minCountsPerGene)

	minSamples=.01*nCells # 1% of cells

	print("minSamples: ", minSamples)

	### adding metadata from seurat schwannoma subcluster object
	schwann_seurat = pd.read_csv(f_metadata, sep='\t', header=0, index_col=0 )
	print(schwann_seurat.head())
	adata.obs['Sample_name'] = schwann_seurat['orig.ident']
	adata.obs['technique'] = schwann_seurat['technique']
	adata.obs['Seurat_clustering'] = schwann_seurat['integrated_snn_res.0.5']
	adata.obs['schwann_cluster_annotation'] = schwann_seurat['schwann_cluster_functionalannot']

	# [10] simply compute the number of genes per cell (computes 'n_genes' column)
	sc.pp.filter_cells(adata, min_genes=0)
	# mito and genes/counts cuts
	mito_genes = adata.var_names.str.startswith('MT-')
	print(mito_genes)

	# for each cell compute fraction of counts in mito genes vs. all genes
	adata.obs['percent_mito'] = np.sum(
	adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
	# add the total counts per cell as observations-annotation to adata
	adata.obs['n_counts'] = adata.X.sum(axis=1).A1

	### carry out filtering steps
	# [16]
	# initial cuts
	sc.pp.filter_cells(adata, min_genes=200 )
	sc.pp.filter_genes(adata, min_cells=3 )

	# [17]
	adata = adata[adata.obs['n_genes'] < 4000, :] # may want to change this
	adata = adata[adata.obs['percent_mito'] < 0.15, :]

	### Finalize the selected filters
	### Update the anndata file, to be used for further processing, clustering, visualization, etc..

	# [21] update the ann data:
	adata.write( f_anndata_path ) # this outputs the results to a loom file

	# [25] create basic row and column attributes for the loom file:
	row_attrs = {
	"Gene": np.array(adata.var_names) ,
	}
	col_attrs = {
	"CellID": np.array(adata.obs_names) ,
	"nGene": np.array( np.sum(adata.X.transpose()>0 , axis=0)).flatten() ,
	"nUMI": np.array( np.sum(adata.X.transpose() , axis=0)).flatten() ,
	"Sample_name": np.array(adata.obs['Sample_name']),
	"technique" : np.array(adata.obs["technique"]),
	"Seurat_clustering" : np.array(adata.obs["Seurat_clustering"]),
	"schwann_cluster_annotation" : np.array(adata.obs["schwann_cluster_annotation"]),
	}

	# Index(['nGene', 'nUMI', 'Sample_name', 'technique', 'Seurat_clustering',
	#        'schwann_cluster_annotation', 'n_genes', 'percent_mito', 'n_counts'],
	#       dtype='object')

	lp.create( f_loom_path_scenic, adata.X.transpose(), row_attrs, col_attrs)


	### further pre-processing of expression data
	# Block [20]: save a copy of the raw data
	adata.raw = adata

	# Total-count normalize (library-size correct) to 10,000 reads/cell
	sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)

	# log transform the data.
	sc.pp.log1p(adata)

	# identify highly variable genes.
	sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
	sc.pl.highly_variable_genes(adata)

	# keep only highly variable genes:
	adata = adata[:, adata.var['highly_variable']]

	# regress out total counts per cell and the percentage of mitochondrial genes expressed
	sc.pp.regress_out(adata, ['n_counts', 'percent_mito'] ) #, n_jobs=args.threads)

	# scale each gene to unit variance, clip values exceeding SD 10.
	sc.pp.scale(adata, max_value=10)

	# update the anndata file:
	adata.write( f_anndata_path )

	# Block [21]: run PCA
	sc.tl.pca(adata, svd_solver='arpack')
	sc.pl.pca_variance_ratio(adata, log=True)
	adata.write( f_anndata_path )

	# Block [22]: neighborhood graph of cells (determine optimal number of PCs here)
	sc.pp.neighbors(adata, n_neighbors=15, n_pcs=30)
	# compute UMAP
	sc.tl.umap(adata)
	# tSNE
	tsne = TSNE( n_jobs=20 )
	adata.obsm['X_tsne'] = tsne.fit_transform( adata.X )
	adata.write( f_anndata_path )

	### Clustering
	# Block [23]: cluster the neighbourhood graph
	sc.tl.louvain(adata,resolution=0.4)

	sc.pl.umap(adata, color=['louvain'] )

	# Block [24]: find marker genes
	sc.tl.rank_genes_groups(adata, 'louvain', method='t-test')
	sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
	# sc.tl.rank_genes_groups(adata, 'louvain', method='logreg')
	# sc.pl.rank_genes_groups(adata, n_genes=25, sharey=False)
	pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(10)
	adata.write( f_anndata_path )
