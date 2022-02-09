
countfn='/storage1/fs1/alberthkim/Active/users/khan.saad/NMF_analysis_iteration2/cNMF_analysis/scwhann_scrna.nmf_out.h5ad'

hvgs = open('/storage1/fs1/alberthkim/Active/users/khan.saad/NMF_analysis_iteration2/cNMF_analysis/schwann_scrna_cnmf/schwann_scrna_cnmf.overdispersed_genes.txt').read().split('\n')


adata = sc.read(countfn)


## Set log-normalized data to the raw attribute of the AnnData object to make it easy to plot expression levels of individual genes.
## This does not log normalize the actual AnnData data matrix

sc.pp.normalize_per_cell(adata, counts_per_cell_after=10**4) ## TPT normalization

output_h5ad='/storage1/fs1/alberthkim/Active/users/khan.saad/NMF_analysis_iteration2/cNMF_analysis/scwhann_scrna.nmf_out_pca_filt_normalized.h5ad'

sc.write(output_h5ad, adata)

adata.raw = sc.pp.log1p(adata.copy(), copy=True)

output_h5ad='/storage1/fs1/alberthkim/Active/users/khan.saad/NMF_analysis_iteration2/cNMF_analysis/scwhann_scrna.nmf_out_pca_filt_normalized_raw.h5ad'

sc.write(output_h5ad, adata)


## Subset out only the high-variance genes

adata = adata[:,hvgs]

## Mean and variance normalize the genes

sc.pp.scale(adata)

## Run PCA

sc.pp.pca(adata)

sc.pl.pca_variance_ratio(adata, log=True)


sc.pl.pca_variance_ratio(adata, log=True, save='pca_plot_schwann.png')

## Construct the nearest neighbor graph for UMAP

sc.pp.neighbors(adata, n_neighbors=50, n_pcs=15)

## Run UMAP

sc.tl.umap(adata)

output_h5ad='/storage1/fs1/alberthkim/Active/users/khan.saad/NMF_analysis_iteration2/cNMF_analysis/scwhann_scrna.nmf_out_pca_filt_normalized_scaled_hvg_umap.h5ad'

sc.write(output_h5ad, adata)



output_h5ad='/storage1/fs1/alberthkim/Active/users/khan.saad/NMF_analysis_iteration2/cNMF_analysis/scwhann_scrna.nmf_out_pca_filt_normalized_scaled_hvg_umap.h5ad'

# sc.read(output_h5ad, adata)

#remove columns from anndata object 
# idx = [0]
# new_col_names = 'n_counts'
# adata.obs[new_col_names] = adata.obs[:,idx]


# output_h5ad='/storage1/fs1/alberthkim/Active/users/khan.saad/NMF_analysis_iteration2/cNMF_analysis/scwhann_scrna.nmf_out_pca_filt_hvgonly.h5ad'

# sc.write(output_h5ad, adata)

## Plot the UMAP with some cannonical marker genes to see that the apparent clustering makes sense

sc.pl.umap(adata, color=['IL7R', 'GZMB', 'CD8A','GNLY',
                         'MS4A1', 'CD14', 'FCGR3A', 'HLA-DRA'],
           use_raw=True, ncols=3,save='umap_genes_schwann.png')

# usage = pd.read_csv('/storage1/fs1/alberthkim/Active/users/khan.saad/NMF_analysis_iteration2/cNMF_analysis/schwann_scrna_cnmf/schwann_scrna_cnmf.usages.k_9.dt_0_10.consensus.txt',sep='\t', index_col=0)
# usage.columns = ['Usage_%s' % i for i in usage.columns]
# usage.head()
usage = pd.read_csv('/storage1/fs1/alberthkim/Active/users/khan.saad/NMF_analysis_iteration2/cNMF_analysis/schwann_scrna_cnmf/schwann_scrna_cnmf.usages.k_10.dt_0_05.consensus.txt',sep='\t', index_col=0)
# usage = pd.read_csv('/storage1/fs1/alberthkim/Active/users/khan.saad/NMF_analysis_iteration2/cNMF_analysis/schwann_scrna_cnmf/schwann_scrna_cnmf.usages.k_9.dt_0_10.consensus.txt',sep='\t', index_col=0)
usage.columns = ['Usage_%s_K10' % i for i in usage.columns]
usage.head()

#/storage1/fs1/alberthkim/Active/users/khan.saad/NMF_analysis_iteration2/cNMF_analysis/schwann_scrna_cnmf/schwann_scrna_cnmf.usages.k_7.dt_0_10.consensus.txt

usage = pd.read_csv('/storage1/fs1/alberthkim/Active/users/khan.saad/NMF_analysis_iteration2/cNMF_analysis/schwann_scrna_cnmf/schwann_scrna_cnmf.usages.k_7.dt_0_10.consensus.txt',sep='\t', index_col=0)
usage.columns = ['Usage_%s' % i for i in usage.columns]
usage.head()

# hvgs = open('/storage1/fs1/alberthkim/Active/users/khan.saad/NMF_analysis_iteration2/cNMF_analysis/schwann_scrna_cnmf/schwann_scrna_cnmf.overdispersed_genes.txt').read().split('\n')

# sc.pp.normalize_per_cell(adata, counts_per_cell_after=10**4) ## TPT normalization

# adata.raw = sc.pp.log1p(adata.copy(), copy=True)

# adata = adata[:,hvgs]

# sc.write(output_h5ad, adata)

# output_h5ad='/storage1/fs1/alberthkim/Active/users/khan.saad/NMF_analysis_iteration2/cNMF_analysis/scwhann_scrna.nmf_out_pca_filt_raw.h5ad'

# sc.write(output_h5ad, adata.raw)


# sc.pp.scale(adata)

# 
# adata.raw = sc.pp.log1p(adata.copy(), copy=True)


usage_norm = usage.div(usage.sum(axis=1), axis=0)

adata.obs = pd.merge(left=adata.obs, right=usage_norm, how='left', left_index=True, right_index=True)

# sc.pl.umap(adata,use_raw=True, color=usage_norm.columns,ncols=3, vmin=0, vmax=1,save='umap_usage_K7_schwann.png')

# 'Usage_2_K10', 'Usage_3_K10',
#        'Usage_4_K10', 'Usage_5_K10', 'Usage_6_K10', 'Usage_7_K10',
#        'Usage_8_K10', 'Usage_9_K10', 'Usage_10_K10'

sc.pl.umap(adata,use_raw=True, color=usage_norm.columns,ncols=3, vmin=0, vmax=1,save='umap_usage_K10_schwann.png')


# gene_scores = pd.read_csv('/storage1/fs1/alberthkim/Active/users/khan.saad/NMF_analysis_iteration2/cNMF_analysis/schwann_scrna_cnmf/schwann_scrna_cnmf.gene_spectra_score.k_7.dt_0_10.txt',sep='\t', index_col=0).T

gene_scores = pd.read_csv('/storage1/fs1/alberthkim/Active/users/khan.saad/NMF_analysis_iteration2/cNMF_analysis/schwann_scrna_cnmf/schwann_scrna_cnmf.gene_spectra_score.k_10.dt_0_05.txt',sep='\t', index_col=0).T


gene_scores.head()

top_genes = []
ngenes = 50
for gep in gene_scores.columns:
	top_genes.append(list(gene_scores.sort_values(by=gep, ascending=False).index[:ngenes]))
top_genes = pd.DataFrame(top_genes, index=gene_scores.columns).T
top_genes

#["ANGPTL4","TCF15","CCL3","CD9","MPZ","FOSB","TMSB10","IFI6","WSB1"]

sc.pl.umap(adata,use_raw=True, color=["IFI6","CD9","CCL3","MPZ","PTN","TCF15","IGFBP7","TMSB10","WSB1","JUNB"],ncols=3,save='umap_topgenes_byprogram_schwann_K10.png')

# sc.pl.umap(adata,use_raw=True, color=["WSB1","MPZ","FOSB","CCL3","TMSB10","HLA-DRB5","CD9"],ncols=3,save='umap_topgenes_byprogram_schwann_K7.png')

top_genes.columns = ['program_%s' % i for i in top_genes.columns]

top_genes.to_csv('top_genes_K10_ngene50.csv', index=False) 


top_genes = []
ngenes = 25
for gep in gene_scores.columns:
	top_genes.append(list(gene_scores.sort_values(by=gep, ascending=False).index[:ngenes]))
top_genes = pd.DataFrame(top_genes, index=gene_scores.columns).T
top_genes

# top_genes.columns = ['program_%s' % i for i in top_genes.columns]

# top_genes.to_csv('top_25_genes_K9.csv', index=False)

top_genes.columns = ['program_%s' % i for i in top_genes.columns]

top_genes.to_csv('top_25_genes_K10.csv', index=False)


# usage_norm.to_csv('usage_scores_by_cell_K7.csv', index=True)

usage_norm.to_csv('usage_scores_by_cell_K10.csv', index=True)


seurat_object_dup <- seurat_object
seurat_object <- AddMetaData(seurat_object, usage_scores_cNMF)

