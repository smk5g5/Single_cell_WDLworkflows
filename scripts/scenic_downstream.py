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

import json
import zlib
import base64
import pickle

#This script will be updated to work with WDL

# Regulon specificity scores (RSS) across predicted cell types
from pyscenic.rss import regulon_specificity_scores
from pyscenic.plotting import plot_rss
import matplotlib.pyplot as plt
from adjustText import adjust_text
import seaborn as sns
from pyscenic.binarization import binarize
from pyscenic.cli.utils import load_signatures
from pyscenic.export import add_scenic_metadata
import matplotlib as mpl



import argparse

def dfToNamedMatrix(df):
	arr_ip = [tuple(i) for i in df.values]
	dtyp = np.dtype(list(zip(df.dtypes.index, df.dtypes)))
	arr = np.array(arr_ip, dtype=dtyp)
	return arr

def palplot(pal, names, colors=None, size=1):
	n = len(pal)
	f, ax = plt.subplots(1, 1, figsize=(n * size, size))
	ax.imshow(np.arange(n).reshape(1, n),
		cmap=mpl.colors.ListedColormap(list(pal)),
		interpolation="nearest", aspect="auto")
	ax.set_xticks(np.arange(n) - .5)
	ax.set_yticks([-.5, .5])
	ax.set_xticklabels([])
	ax.set_yticklabels([])
	colors = n * ['k'] if colors is None else colors
	for idx, (name, color) in enumerate(zip(names, colors)):
		ax.text(0.0+idx, 0.0, name, color=color, horizontalalignment='center', verticalalignment='center')
	return f


def rss_regulon_scores(auc_mat,cellannot_dat,output_figure_file):
	rss_cellType = regulon_specificity_scores( auc_mat, cellannot_dat)
	cats = sorted(list(set(cellannot_dat)))

	fig = plt.figure(figsize=(15, 8))
	for c,num in zip(cats, range(1,len(cats)+1)):
		x=rss_cellType.T[c]
		if len(cats) <=10:
			ax = fig.add_subplot(2,5,num)
		else:
			ax = fig.add_subplot(3,4,num)
		plot_rss(rss_cellType, c, top_n=5, max_n=None, ax=ax)
		ax.set_ylim( x.min()-(x.max()-x.min())*0.05 , x.max()+(x.max()-x.min())*0.05 )
		for t in ax.texts:
			t.set_fontsize(12)
		ax.set_ylabel('')
		ax.set_xlabel('')
		adjust_text(ax.texts, autoalign='xy', ha='right', va='bottom', arrowprops=dict(arrowstyle='-',color='lightgrey'), precision=0.001 )

	fig.text(0.5, 0.0, 'Regulon', ha='center', va='center', size='x-large')
	fig.text(0.00, 0.5, 'Regulon specificity score (RSS)', ha='center', va='center', rotation='vertical', size='x-large')
	plt.tight_layout()
	plt.rcParams.update({
	'figure.autolayout': True,
	    'figure.titlesize': 'large' ,
	    'axes.labelsize': 'medium',
	    'axes.titlesize':'large',
	    'xtick.labelsize':'medium',
	    'ytick.labelsize':'medium'
	    })
	# plt.savefig("Schwann_subclus_louvain_clustering-RSS-top5.pdf", dpi=600, bbox_inches = "tight")
	# plt.tight_layout()

	plt.savefig(output_figure_file, bbox_inches = 'tight')
	return(rss_cellType)

def plot_heatmap(auc_mat,cellannot_dat,output_meta,tech,auc_zscore):
	cats = sorted(list(set(cellannot_dat)))
	output_figure_file = str(tech+'_'+output_meta+'-RSS-top5.pdf')
	rss_cellType = rss_regulon_scores(auc_mat,cellannot_dat,output_figure_file)
	topreg = []
	for i,c in enumerate(cats):
		topreg.extend(list(rss_cellType.T[c].sort_values(ascending=False)[:5].index))
	topreg = list(set(topreg))
	print(topreg)

	colors = sns.color_palette('bright',n_colors=len(cats) )
	colorsd = dict( zip( cats, colors ))
	colormap = [ colorsd[x] for x in cellannot_dat]

	output_heatmap_leg =  str(tech+'_'+output_meta+'_'+'-heatmap-legend-top5.pdf')

	sns.set()
	sns.set(font_scale=0.8)
	fig = palplot( colors, cats, size=1.0)
	plt.savefig(output_heatmap_leg, dpi=300, bbox_inches = "tight")

	output_heatmap =  str(tech+'_'+output_meta+'_'+'-heatmap-top5.pdf')

	sns.set(font_scale=1.2)
	g = sns.clustermap(auc_mtx_Z[topreg], annot=False,  square=False,  linecolor='gray',
	yticklabels=False, xticklabels=True, vmin=-2, vmax=6, row_colors=colormap,
	cmap="YlGnBu", figsize=(21,16) )
	g.cax.set_visible(True)
	g.ax_heatmap.set_ylabel('')
	g.ax_heatmap.set_xlabel('')
	plt.savefig(output_heatmap,bbox_inches = "tight")
	heatmap_dict = dict()
	heatmap_dict['topreg'] = topreg
	heatmap_dict['rss_cellType'] = rss_cellType
	return(heatmap_dict)




if __name__ == "__main__":
	parser=argparse.ArgumentParser(description="""
	This script will run Stream bifurcation analysis
	pass Counts file metadata file output dir and sample name
	""")
	parser.add_argument("scenic_loom", help="scenic loom file")
	parser.add_argument("scanpy_loom", help="scenic loom file")
	parser.add_argument("anndata_input", help="anndata input file")
	parser.add_argument("outdir", help="output directory")
	parser.add_argument("metaname", help="output_file_meta name")

	args = parser.parse_args()

	scenic_loom=args.scenic_loom
	scanpy_loom = args.scanpy_loom
	anndata_input = args.anndata_input
	OUTDIR=args.outdir
	metaname=args.metaname

	wdir = OUTDIR
	os.chdir( wdir )

	f_pyscenic_output = str(wdir+'/'+scenic_loom)
	f_final_loom = str(wdir+'/'+ metaname + '.scenic.scope.loom')
	f_loom_path_scenic = str(wdir+'/'+ scanpy_loom)
	f_anndata_file = str(wdir+'/'+ anndata_input)

	print('f_pyscenic_output',f_pyscenic_output)
	print('f_final_loom',f_final_loom)
	print('f_loom_path_scenic',f_loom_path_scenic)
	print('f_anndata_file',f_anndata_file)

	# collect SCENIC AUCell output
	lf = lp.connect( f_pyscenic_output, mode='r+', validate=False )
	auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)

	import umap
	runUmap = umap.UMAP(n_neighbors=10, min_dist=0.4,metric='correlation').fit_transform
	dr_umap = runUmap(auc_mtx)
	pd.DataFrame(dr_umap, columns=['X', 'Y'], index=auc_mtx.index).to_csv( "scenic_umap.txt", sep='\t') # import umap

	# tSNE
	tsne = TSNE(n_jobs=20)
	dr_tsne = tsne.fit_transform(auc_mtx)
	pd.DataFrame(dr_tsne, columns=['X', 'Y'], index=auc_mtx.index).to_csv( "scenic_tsne.txt", sep='\t')

	meta = json.loads(zlib.decompress(base64.b64decode( lf.attrs.MetaData )))
	#exprMat = pd.DataFrame( lf[:,:], index=lf.ra.Gene, columns=lf.ca.CellID)
#	auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)
	regulons = lf.ra.Regulons
	dr_umap = pd.read_csv( 'scenic_umap.txt', sep='\t', header=0, index_col=0 )
	dr_tsne = pd.read_csv( 'scenic_tsne.txt', sep='\t', header=0, index_col=0 )

	auc_mtx.columns = auc_mtx.columns.str.replace('\(','_(')
	regulons.dtype.names = tuple( [ x.replace("(","_(") for x in regulons.dtype.names ] )

	# regulon thresholds
	rt = meta['regulonThresholds']
	for i,x in enumerate(rt):
	    tmp = x.get('regulon').replace("(","_(")
	    x.update( {'regulon': tmp} )

	adata = sc.read(f_anndata_file)

	print('length of lf.ca.CellID',len(lf.ca.CellID))

	print('shape of adata',adata.shape)

	print('adata.obs.index',adata.obs.index)

	Embeddings_X = pd.DataFrame( index=lf.ca.CellID )
	Embeddings_X = pd.concat( [
	pd.DataFrame(adata.obsm['X_umap'],index=adata.obs.index)[0] ,
	pd.DataFrame(adata.obsm['X_pca'],index=adata.obs.index)[0] ,
	dr_tsne['X'] ,
	dr_umap['X']
	], sort=False, axis=1, join='outer' )
	Embeddings_X.columns = ['1','2','3','4']

	Embeddings_Y = pd.DataFrame( index=lf.ca.CellID )
	Embeddings_Y = pd.concat( [
	pd.DataFrame(adata.obsm['X_umap'],index=adata.obs.index)[1] ,
	pd.DataFrame(adata.obsm['X_pca'],index=adata.obs.index)[1] ,
	dr_tsne['Y'] ,
	dr_umap['Y']
	], sort=False, axis=1, join='outer' )
	Embeddings_Y.columns = ['1','2','3','4']

	adata.raw = adata
	adata = adata[list(Embeddings_X.index)] #subset anndata object so it contains the same cells as f_pyscenic_output file
	
	# concatenate embeddings:
	tsneDF = pd.DataFrame(adata.obsm['X_tsne'], columns=['_X', '_Y'])


	metaJson = {}
	
	metaJson['embeddings'] = [
	    {
	        "id": -1,
	        "name": f"Scanpy t-SNE (highly variable genes)"
	    },
	    {
	        "id": 1,
	        "name": f"Scanpy UMAP  (highly variable genes)"
	    },
	    {
	        "id": 2,
	        "name": "Scanpy PC1/PC2"
	    },
	    {
	        "id": 3,
	        "name": "SCENIC AUC t-SNE"
	    },
	    {
	        "id": 4,
	        "name": "SCENIC AUC UMAP"
	    },
	]
	
	metaJson["clusterings"] = [{
	    "id": 0,
	    "group": "Scanpy",
	    "name": "Scanpy louvain default resolution",
	    "clusters": [],
	}]
	
	metaJson["metrics"] = [
	    {
	            "name": "nUMI"
	    }, {
	            "name": "nGene"
	    }, {
	            "name": "Percent_mito"
	    }
	]
	
	metaJson["annotations"] = [
	    {
	        "name": "Louvain_clusters_Scanpy",
	        "values": list(set( adata.obs['louvain'].astype(np.str) ))
	    },
	    #{
	    #    "name": "Genotype",
	    #    "values": list(set(adata.obs['Genotype'].values))
	    #},
	    #{
	    #    "name": "Timepoint",
	    #    "values": list(set(adata.obs['Timepoint'].values))
	    #},
	    #{
	    #    "name": "Sample",
	    #    "values": list(set(adata.obs['Sample'].values))
	    #}
	]
	
	# SCENIC regulon thresholds:
	metaJson["regulonThresholds"] = rt

	for i in range(max(set([int(x) for x in adata.obs['louvain']])) + 1):
		clustDict = {}
		clustDict['id'] = i
		clustDict['description'] = f'Unannotated Cluster {i + 1}'
		metaJson['clusterings'][0]['clusters'].append(clustDict)

	clusterings = pd.DataFrame()
	clusterings["0"] = adata.obs['louvain'].values.astype(np.int64)


	col_attrs = {
	"CellID": np.array(adata.obs.index),
	"nUMI": np.array(adata.obs['n_counts'].values),
	"nGene": np.array(adata.obs['n_genes'].values),
	"Louvain_clusters_Scanpy": np.array( adata.obs['louvain'].values ),
	#"Genotype": np.array(adata.obs['Genotype'].values),
	#"Timepoint": np.array(adata.obs['Timepoint'].values),
	#"Sample": np.array(adata.obs['Sample'].values),
	"Percent_mito": np.array(adata.obs['percent_mito'].values),
	'Sample_name': np.array(adata.obs['Sample_name'].values),
#	'technique': np.array(adata.obs['technique'].values),
#	'Seurat_clustering': np.array(adata.obs['Seurat_clustering'].values),
	'cluster_annotation': np.array(adata.obs['cluster_annotation'].values),
	"Embedding": dfToNamedMatrix(tsneDF),
	"Embeddings_X": dfToNamedMatrix(Embeddings_X),
	"Embeddings_Y": dfToNamedMatrix(Embeddings_Y),
	"RegulonsAUC": dfToNamedMatrix(auc_mtx),
	"Clusterings": dfToNamedMatrix(clusterings),
	"ClusterID": np.array(adata.obs['louvain'].values)
	}


	row_attrs = {
	"Gene": lf.ra.Gene,
	"Regulons": regulons,
	}

	attrs = {
	"title": "sampleTitle",
	"MetaData": json.dumps(metaJson),
	"Genome": 'hg38',
	"SCopeTreeL1": "",
	"SCopeTreeL2": "",
	"SCopeTreeL3": ""
	}

	attrs['MetaData'] = base64.b64encode(zlib.compress(json.dumps(metaJson).encode('ascii'))).decode('ascii')

	lp.create(
	filename = f_final_loom ,
	layers=lf[:,:],
	row_attrs=row_attrs,
	col_attrs=col_attrs,
	file_attrs=attrs
	)
	lf.close()

	lf = lp.connect( f_final_loom, mode='r', validate=False )
	meta = json.loads(zlib.decompress(base64.b64decode( lf.attrs.MetaData )))
	exprMat = pd.DataFrame( lf[:,:], index=lf.ra.Gene, columns=lf.ca.CellID).T
	auc_mtx = pd.DataFrame( lf.ca.RegulonsAUC, index=lf.ca.CellID)

	# create a dictionary of regulons:
	regulons = {}
	for i,r in pd.DataFrame(lf.ra.Regulons,index=lf.ra.Gene).iteritems():
		regulons[i] =  list(r[r==1].index.values)


	# ['CellID', 'ClusterID', 'Clusterings', 'Embedding', 'Embeddings_X', 'Embeddings_Y', 'Louvain_clusters_Scanpy',
	# 'Percent_mito', 'RegulonsAUC',
	# 'Sample_name', 'Seurat_clustering', 'nGene', 'nUMI', 'schwann_cluster_annotation', 'technique']
	# cell annotations from the loom column attributes:
	cellAnnot = pd.concat(
	[
	    pd.DataFrame( lf.ca.ClusterID, index=lf.ca.CellID ),
	    pd.DataFrame( lf.ca.Louvain_clusters_Scanpy, index=lf.ca.CellID ),
	    pd.DataFrame( lf.ca.Percent_mito, index=lf.ca.CellID ),
	    pd.DataFrame( lf.ca.nGene, index=lf.ca.CellID ),
	    pd.DataFrame( lf.ca.nUMI, index=lf.ca.CellID ),
	    pd.DataFrame( lf.ca.Sample_name, index=lf.ca.CellID ),
#	    pd.DataFrame( lf.ca.Seurat_clustering, index=lf.ca.CellID ),
#	    pd.DataFrame( lf.ca.technique, index=lf.ca.CellID ),
	    pd.DataFrame( lf.ca.cluster_annotation, index=lf.ca.CellID ),
	],
	axis=1
	)

	cellAnnot.columns = [
	'ClusterID',
	'Louvain_clusters_Scanpy',
	'Percent_mito',
	'nGene',
	'nUMI',
	'Sample_name',
	'cluster_annotation']
#	'Seurat_clustering',
#	'technique',
#	'cluster_annotation'

	 # capture embeddings:
	dr = [
	pd.DataFrame( lf.ca.Embedding, index=lf.ca.CellID )
	]
	dr_names = [
	meta['embeddings'][0]['name'].replace(" ","_")
	]

	# add other embeddings
	drx = pd.DataFrame( lf.ca.Embeddings_X, index=lf.ca.CellID )
	dry = pd.DataFrame( lf.ca.Embeddings_Y, index=lf.ca.CellID )


	for i in range( len(drx.columns) ):
		dr.append( pd.concat( [ drx.iloc[:,i], dry.iloc[:,i] ], sort=False, axis=1, join='outer' ))
		dr_names.append( meta['embeddings'][i+1]['name'].replace(" ","_").replace('/','-') )

	# rename columns:
	for i,x in enumerate( dr ):
		x.columns = ['X','Y']

	lf.close()

	for i,x in enumerate( dr ):
		adata.obsm[ 'X_'+dr_names[i] ] = x.to_numpy()


	sc._utils.sanitize_anndata( adata )

	sig = load_signatures('reg.csv')
	adata = add_scenic_metadata(adata, auc_mtx, sig)

	from pyscenic.utils import load_motifs
	import operator as op
	from IPython.display import HTML, display

	BASE_URL = "http://motifcollections.aertslab.org/v9/logos/"
	COLUMN_NAME_LOGO = "MotifLogo"
	COLUMN_NAME_MOTIF_ID = "MotifID"
	COLUMN_NAME_TARGETS = "TargetGenes"

	from pyscenic.utils import load_motifs
	import operator as op
	from IPython.display import HTML, display

	BASE_URL = "http://motifcollections.aertslab.org/v9/logos/"
	COLUMN_NAME_LOGO = "MotifLogo"
	COLUMN_NAME_MOTIF_ID = "MotifID"
	COLUMN_NAME_TARGETS = "TargetGenes"

	def display_logos(df: pd.DataFrame, top_target_genes: int = 3, base_url: str = BASE_URL):
		"""
		:param df:
		:param base_url:
		"""
		# Make sure the original dataframe is not altered.
		df = df.copy()

		# Add column with URLs to sequence logo.
		def create_url(motif_id):
			return '<img src="{}{}.png" style="max-height:124px;"></img>'.format(base_url, motif_id)
		df[("Enrichment", COLUMN_NAME_LOGO)] = list(map(create_url, df.index.get_level_values(COLUMN_NAME_MOTIF_ID)))

		# Truncate TargetGenes.
		def truncate(col_val):
			return sorted(col_val, key=op.itemgetter(1))[:top_target_genes]
		df[("Enrichment", COLUMN_NAME_TARGETS)] = list(map(truncate, df[("Enrichment", COLUMN_NAME_TARGETS)]))

		MAX_COL_WIDTH = pd.get_option('display.max_colwidth')
		pd.set_option('display.max_colwidth', 200)
		a = HTML(df.head().to_html(escape=False))
		print(a)
		pd.set_option('display.max_colwidth', MAX_COL_WIDTH)
		with open('motif_output.html', 'w') as f:
			f.write(a)

	df_motifs = load_motifs('reg.csv')

	selected_motifs = ['PAX5','TCF3','EBF1']
	df_motifs_sel = df_motifs.iloc[ [ True if x in selected_motifs else False for x in df_motifs.index.get_level_values('TF') ] ,:]

	sc.set_figure_params(frameon=False, dpi=600, fontsize=10, dpi_save=600)

	sc.pl.scatter( adata, basis='Scanpy_UMAP__(highly_variable_genes)',
	color=['louvain','Sample_name','cluster_annotation'],
	title=['HVG - UMAP (Louvain clusters)','HVG - UMAP (Sample name)'],
	alpha=0.8,
	save='_Louvain-celltype.pdf'
	)

	sc.pl.scatter( adata, basis='SCENIC_AUC_UMAP',
	color=['louvain','Sample_name','cluster_annotation'],
	title=['SCENIC - UMAP (Louvain clusters)','SCENIC - UMAP (Sample name)','SCENIC - UMAP (cluster_annotation)'],
	alpha=0.8,
	save='_Louvain-celltype.pdf'
	)

	#### Generate a Z-score for each regulon to enable comparison between regulons
	auc_mtx_Z = pd.DataFrame( index=auc_mtx.index )
	for col in list(auc_mtx.columns):
		auc_mtx_Z[col] = (auc_mtx[col] - auc_mtx[col].mean()) / auc_mtx[col].std(ddof=0)

	output_dict = dict()

	output_dict['Louvain_clusters_Scanpy'] = plot_heatmap(auc_mat=auc_mtx,cellannot_dat=cellAnnot['Louvain_clusters_Scanpy'],output_meta='Louvain_clusters_Scanpy',tech=metaname,auc_zscore=auc_mtx_Z)
	output_dict['cluster_annotation'] = plot_heatmap(auc_mat=auc_mtx,cellannot_dat=cellAnnot['cluster_annotation'],output_meta='cluster_annotation',tech=metaname,auc_zscore=auc_mtx_Z)
	output_dict['Sample_name'] = plot_heatmap(auc_mat=auc_mtx,cellannot_dat=cellAnnot['Sample_name'],output_meta='Sample_name',tech=metaname,auc_zscore=auc_mtx_Z)

	pickle_output = str(metaname+'_'+'dictobj.pkl')
	filehandler = open(pickle_output,"wb")
	pickle.dump(output_dict,filehandler)
	filehandler.close()
