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
'metalabels': np.array(adata.obs['metalabels'].values),
#	'Seurat_clustering': np.array(adata.obs['Seurat_clustering'].values),
'label': np.array(adata.obs['label'].values),
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


regulons = {}
for i,r in pd.DataFrame(lf.ra.Regulons,index=lf.ra.Gene).iteritems():
	regulons[i] =  list(r[r==1].index.values)


cellAnnot = pd.concat(
[
    pd.DataFrame( lf.ca.ClusterID, index=lf.ca.CellID ),
    pd.DataFrame( lf.ca.Louvain_clusters_Scanpy, index=lf.ca.CellID ),
    pd.DataFrame( lf.ca.Percent_mito, index=lf.ca.CellID ),
    pd.DataFrame( lf.ca.nGene, index=lf.ca.CellID ),
    pd.DataFrame( lf.ca.nUMI, index=lf.ca.CellID ),
    pd.DataFrame( lf.ca.Sample_name, index=lf.ca.CellID ),
#	    pd.DataFrame( lf.ca.Seurat_clustering, index=lf.ca.CellID ),
	pd.DataFrame( lf.ca.metalabels, index=lf.ca.CellID ),
    pd.DataFrame( lf.ca.label, index=lf.ca.CellID ),
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
'metalabels',
'label']

dr = [
pd.DataFrame( lf.ca.Embedding, index=lf.ca.CellID )
]
dr_names = [
meta['embeddings'][0]['name'].replace(" ","_")
]

drx = pd.DataFrame( lf.ca.Embeddings_X, index=lf.ca.CellID )
dry = pd.DataFrame( lf.ca.Embeddings_Y, index=lf.ca.CellID )



for i in range( len(drx.columns) ):
	dr.append( pd.concat( [ drx.iloc[:,i], dry.iloc[:,i] ], sort=False, axis=1, join='outer' ))
	dr_names.append( meta['embeddings'][i+1]['name'].replace(" ","_").replace('/','-') )


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

auc_mtx_Z = pd.DataFrame( index=auc_mtx.index )
for col in list(auc_mtx.columns):
	auc_mtx_Z[col] = (auc_mtx[col] - auc_mtx[col].mean()) / auc_mtx[col].std(ddof=0)

output_dict['Louvain_clusters_Scanpy'] = plot_heatmap(auc_mat=auc_mtx,cellannot_dat=cellAnnot['Louvain_clusters_Scanpy'],output_meta='Louvain_clusters_Scanpy',tech=metaname,auc_zscore=auc_mtx_Z)



label

metalabels

st.plot_flat_tree(adata,color=['label'],
        dist_scale=1,show_graph=True,show_text=True,
        save_fig=True,fig_size=(8,8),fig_name='flat_tree_label.png',
        vmin=0,vmax=1)

st.detect_transition_markers

use_precomputed=False


st.detect_transition_markers(adata,marker_list=[x.name for x in adata.uns['scaled_marker_expr']],
                             cutoff_spearman=0.4,cutoff_logfc=0.25,use_precomputed=True, root='S6',n_jobs=8)


/storage1/fs1/alberthkim/Active/users/khan.saad/STREAM_analysis/Stream_pickle/stream_result_bifurcation_schwann_nd_fibroblasts_SCH9_adata.pkl

/opt/conda/bin/python


import stream as st
print(st.__version__)
import os
os.getcwd()
import argparse

st.set_figure_params(dpi=300,style='white',figsize=[5.4,4.8],
rc={'image.cmap': 'viridis'})
adata=st.read(file_name=str(PICKLEINP))
st.set_workdir(adata,workdir=str(OUTDIR))
pseudo_times = list(filter(lambda x:'_pseudotime' in x, list(adata.obs.columns)))
st.plot_visualization_2D(adata,n_neighbors=50,color=['label'],
save_fig=True,fig_size=(10,10),fig_name='plot_visualization_epg_2D_label.png')
st.plot_visualization_2D(adata,n_neighbors=50,color=['Sample'],
save_fig=True,fig_size=(10,10),fig_name='plot_visualization_epg_2D_sample.png')
st.plot_visualization_2D(adata,n_neighbors=50,color=['metalabels'],
save_fig=True,fig_size=(10,10),fig_name='plot_visualization_epg_2D_clusters.png')
st.plot_visualization_2D(adata,n_neighbors=50,color=['branch_id_alias'],
save_fig=True,fig_size=(10,10),fig_name='plot_visualization_epg_2D_branchid.png')


pseudo_times = list(filter(lambda x:'_pseudotime' in x, list(adata_allsamp.obs.columns)))

for i in pseudo_times:
	print(i)
	myfig='flat_tree_changedscale_'+str(i)+'.html'
	myfig2='extended_elastic_principal_graph_changedscale_'+str(i)+'.html'
	flat_tree=st.plot_flat_tree(adata_allsamp,color=[i],dist_scale=2,show_graph=True,show_text=True,plotly=True,fig_size=(10,10),vmin=min(adata_allsamp.obs[i]),vmax=max(adata_allsamp.obs[i]))
	dim_red=st.plot_dimension_reduction(adata_allsamp,color=[i],n_components=3,show_graph=True,show_text=True,plotly=True,fig_size=(10,10),vmin=min(adata_allsamp.obs[i]),vmax=max(adata_allsamp.obs[i]))
	flat_tree.write_html(myfig)
	dim_red.write_html(myfig2)




st.set_figure_params(dpi=600,style='white',figsize=[10,10],palette=sns.color_palette('Spectral').reverse())



