import anndata2ri
from rpy2.robjects import r
anndata2ri.activate()
import scib
import pandas as pd
import scanpy as sc
from datetime import date
import anndata
import random

today = date.today()

sc.set_figure_params(dpi=600)
sc.settings.verbosity = 1
sc.set_figure_params(figsize=(10, 10))
sc.logging.print_header()

# from matplotlib import pyplot as plt

# with plt.rc_context():  # Use this to set figure params like size and dpi
#     sc.pl.plotting_function(..., show=False)
#     plt.savefig("path/to/file.extension", bbox_inches="tight")

# with rc_context({'figure.figsize': (4, 4)}):
#     sc.pl.umap(pbmc, color='CD79A')

def phase_without_regressing_by_cellcycle(anndata,phase_pal=phase_pal,s_genes=s_genes,g2m_genes=g2m_genes):
	#without regressing out cell cycle
	sc.pp.normalize_total(anndata, target_sum=1e6)
	sc.pp.log1p(anndata)
	sc.pp.scale(anndata)
	sc.tl.score_genes_cell_cycle(anndata, s_genes=s_genes, g2m_genes=g2m_genes)
	#calculate pca on anndata
	sc.tl.pca(anndata)
	fig_out = 'Before_cellcycle_regression_removecellswithzerocounts'+'.png'
	sc.pl.pca_scatter(anndata, color='phase',palette=phase_pal,save=fig_out)

anndata  = scib.preprocessing.read_seurat('/storage1/fs1/allegra.petti/Active/GBM/CellState/GBM.9/Redo_PCA_clustering/GBM_cellstate_selgenes.220516.rds')

anndata.raw = anndata  # at the point during preprocessing at which you wish store a copy for visualization and differential testing
sc.pp.filter_cells(anndata, min_counts=1)

Organism='human'

if(str(Organism))=='human':
	cell_cycle_df = pd.read_csv("/storage1/fs1/allegra.petti/Active/10xGenomics/key.gene.lists/CellCycleTirosh.txt",sep="\t",header=0)
	s_genes = list(cell_cycle_df['Gene.Symbol'][cell_cycle_df['List']=='G1/S'])
	g2m_genes = list(cell_cycle_df['Gene.Symbol'][cell_cycle_df['List']=='G2/M'])
	cell_cycle_genes = list(cell_cycle_df['Gene.Symbol'])
	cell_cycle_genes = [x for x in cell_cycle_genes if x in anndata.var_names]
	# cell_cycle_genes = [x.strip() for x in open('/data/regev_lab_cell_cycle_genes.txt')]
else:
	cell_cycle_df = pd.read_csv("/storage1/fs1/allegra.petti/Active/10xGenomics/key.gene.lists/CellCycleTirosh_mouse.txt",sep="\t",header=None)
	cell_cycle_df =['List', 'Gene.Symbol']
	s_genes = list(cell_cycle_df['Gene.Symbol'][cell_cycle_df['List']=='G1/S'])
	g2m_genes = list(cell_cycle_df['Gene.Symbol'][cell_cycle_df['List']=='G2/M'])

phase_pal = {
	'G1':"#4477AA",
	'G2M':"#DDCC77",
	'S' : "#CC6677"
}

phase_without_regressing_by_cellcycle(anndata)

anndata  = scib.preprocessing.read_seurat('/storage1/fs1/allegra.petti/Active/GBM/CellState/GBM.9/Redo_PCA_clustering/GBM_cellstate_selgenes.220516.rds')
sc.pp.filter_cells(anndata, min_counts=1)

phase_without_regressing_by_cellcycle(anndata)
#with regressing cell cycle
anndata  = scib.preprocessing.read_seurat('/storage1/fs1/allegra.petti/Active/GBM/CellState/GBM.9/Redo_PCA_clustering/GBM_cellstate_selgenes.220516.rds')
anndata.raw = anndata  # at the point during preprocessing at which you wish store a copy for visualization and differential testing
sc.pp.filter_cells(anndata, min_counts=1)

sc.tl.score_genes_cell_cycle(anndata, s_genes=s_genes, g2m_genes=g2m_genes)
# Regress out effects of cellcycle
sc.tl.pca(anndata)

anndata.var['mt'] = anndata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(anndata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)


fig_out = '_multipanel_QC'+'.png'

sc.pl.violin(anndata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'],
             jitter=0.4, multi_panel=True,save=fig_out)

fig_out = '_topgenes'+'.png'
#top 20 genes
sc.pl.highest_expr_genes(anndata, n_top=20,save=fig_out)

fig_out = '_total_counts_vs_pct_counts_mt'+'.png'

sc.pl.scatter(anndata, x='total_counts', y='pct_counts_mt',save=fig_out)

fig_out = '_total_counts_vs_n_genes_by_counts'+'.png'

sc.pl.scatter(anndata, x='total_counts', y='n_genes_by_counts',save=fig_out)

sc.pp.normalize_total(anndata, target_sum=1e6)
sc.pp.log1p(anndata)
sc.pp.highly_variable_genes(anndata, min_mean=0.0125, max_mean=3, min_disp=0.5)

fig_out = '_HVGs'+'.png'
sc.pl.highly_variable_genes(anndata,save=fig_out)

sc.pp.regress_out(anndata, ['S_score', 'G2M_score'],n_jobs=10)
sc.pp.scale(anndata)
fig_out = '_after_cellcycle_regression'+'.png'
sc.pl.pca_scatter(anndata, color='phase',palette=phase_pal,save=fig_out)

sc.tl.pca(anndata, svd_solver='arpack')

fig_out = '_after_cellcycle_regression_Principal_component_variance'+'.png'

sc.pl.pca_variance_ratio(anndata, log=True,save=fig_out)

sc.pp.neighbors(anndata, n_neighbors=10, n_pcs=50)
sc.tl.leiden(anndata)

sc.tl.paga(anndata)
sc.pl.paga(anndata, plot=False)  # remove `plot=False` if you want to see the coarse-grained graph
sc.tl.umap(anndata, init_pos='paga')

sc.tl.umap(anndata)

fig_out = '_after_cellcycle_regression_selgenes_UMAP'+'.png'

sc.pl.umap(anndata, color=['PTPRC', 'NKG7', 'TREM2'],save=fig_out)

# As we set the .raw attribute of adata, the previous plots 
# showed the “raw” (normalized, logarithmized, but uncorrected)
#  gene expression. You can also plot the scaled and corrected 
#  gene expression by explicitly stating that you don’t want to use .raw.

fig_out = '_after_cellcycle_regression_selgenes_noraw'+'.png'


sc.pl.umap(anndata, color=['PTPRC', 'NKG7', 'TREM2'], use_raw=False,save=fig_out)


fig_out = '_after_cellcycle_regression_selgenes_clus'+'.png'


sc.pl.umap(anndata, color=['leiden','PTPRC', 'NKG7', 'TREM2'],save=fig_out)

fig_out = '_QC_umapplots'+'.png'

sc.pl.umap(anndata, color=['n_genes_by_counts','total_counts', 'pct_counts_mt'],save=fig_out)


# phase_pal = {
# 	'G1':"#4477AA",
# 	'G2M':"#DDCC77",
# 	'S' : "#CC6677"
# }

fig_out = '_cellcycle_umapplots'+'.png'

sc.pl.umap(anndata,color='phase',palette=phase_pal,save=fig_out)


sample_palette = {
	'B150':"#88CCEE",
	'B152':"#CC6677",
	'B178' : "#DDCC77",
	'B183' : "#117733",
	'B185' : "#332288",
	'B186' : "#AA4499",
	'ST073021' : "#44AA99" ,
	'WU1225_core' : "#999933",
	'WU1225_edge' : "#882255",
	'WU1226_core_CD45neg' : "#661100",
	'WU1226_edge_CD45neg' : "#888888",
}

fig_out = '_sample_umapplots'+today.strftime("%m%d%y")+'.png'

sc.pl.umap(anndata,color='orig.ident',palette=sample_palette,save=fig_out,s=50)

fig_out = '_SingleR_scsorter_umapplots'+'.png'

sc.pl.umap(anndata,color='Cell_Types_SingleR_scsorter_v2',save=fig_out)

# df.colName.value_counts()

        # B150                B152                B178                B183
        #   "#88CCEE"           "#CC6677"           "#DDCC77"           "#117733"
        #        B185                B186            ST073021         WU1225_core
        #   "#332288"           "#AA4499"           "#44AA99"           "#999933"
        # WU1225_edge WU1226_core_CD45neg WU1226_edge_CD45neg
        #   "#882255"           "#661100"           "#888888"

results_file = 'Scanpy_end2end_analysis/GBM.CellState.' + today.strftime("%m%d%y") + '.h5ad'  # the file that will store the analysis results


anndata.write(results_file)


anndata = anndata.read_h5ad('/storage1/fs1/allegra.petti/Active/GBM/CellState/GBM.9/Scanpy_end2end_analysis/GBM.CellState.051722.h5ad')

anndata.obs['Cell_Types_SingleR_scsorter_v2'].value_counts()


celltype_colorpalette = {
	'nef.Malignant.NPC1':"#fb298a",
	'nef.Malignant.MES1':"#2b9b00",
	'nef.Malignant.AC' : "#f846ba",
	'nef.Malignant.OPC' : "#01ca75",
	'br_imm_atl.Mg-TAM' : "#ff84f5",
	'br_imm_atl.Mo-TAM' : "#006a16",
	'nef.Malignant.MES2' : "#fc3557",
	'nef.Malignant.NPC2' : "#00cb90",
	'br_imm_atl.T cells' : "#ff5a76",
	'br_imm_atl.prol. TAM' : "#00c7d9",
	'stb.NPC' : "#f1552b",
	'br_imm_atl.Monocytes' : "#008db1",
	'No.Prediction' : "#e0e0e0",
	'nef.Malignant.Unknown' : "#c0aaff",
	'br_imm_atl.DC' : "#d6c94a",
	'br_imm_atl.NK cells' : "#324f97",
	'stb.InN1b' : "#793882",
	'stb.ExN1_4' : "#9dd48d",
	'stb.Micro' : "#a31925",
	'stb.Astro1' : "#018b6a",
	'stb.Perc' : "#ff8670",
	'br_imm_atl.B cells' : "#185f2a",
	'stb.OPC1' : "#8b3358",
	'stb.Astro4' : "#5c6100",
	'stb.Endo' : "#ffaa5d",
	'stb.Olig3' : "#8f360e",
	'stb.InN5' : "#ffb377",
}

fig_out = '_SingleR_scsorter_umapplots.'+today.strftime("%m%d%y")+'.png'

# s=50, frameon=False, ncols=4,
# ,palette=sample_palette,save=fig_out

sc.set_figure_params(dpi=300)
sc.settings.verbosity = 1
sc.set_figure_params(figsize=(20, 20))
sc.set_figure_params(scanpy=True, fontsize=14)



sc.pl.umap(anndata,color='Cell_Types_SingleR_scsorter_v2',palette=celltype_colorpalette,s=50,save=fig_out)



sc.pl.umap(anndata,color='orig.ident',palette=sample_palette,save=fig_out,s=50)


fig_out = '_SingleR_scsorter_umapplots_grouped'+today.strftime("%m%d%y")+'.png'


def cluster_small_multiples(adata, clust_key,figout,palette,size=50, frameon=False, legend_loc=None,**kwargs):
	tmp = adata.copy()
	for i,clust in enumerate(adata.obs[clust_key].cat.categories):
		tmp.obs[clust] = adata.obs[clust_key].isin([clust]).astype('category')
		tmp.uns[clust+'_colors'] = ['#d3d3d3', adata.uns[clust_key+'_colors'][i]]
		ncol_arg = min(len(adata.obs[clust_key].cat.categories.tolist()),3)
	sc.pl.umap(tmp, groups=tmp.obs[clust].cat.categories[1:].values,color=adata.obs[clust_key].cat.categories.tolist(),ncols=ncol_arg,save=fig_out,size=size, frameon=frameon, legend_loc=legend_loc, **kwargs)

# sc.pl.umap(anndata,color='Cell_Types_SingleR_scsorter_v2',palette=celltype_colorpalette,s=50,save=fig_out,frameon=False, ncols=4,groups=list(celltype_colorpalette.keys()))

fig_out = '_sample_umapplots_grouped'+today.strftime("%m%d%y")+'.png'

# sc.pl.umap(anndata, groups=[[c] for c in anndata.obs['Cell_Types_SingleR_scsorter_v2'].cat.categories],save=fig_out,color='Cell_Types_SingleR_scsorter_v2',palette=celltype_colorpalette, ncols=4)#, simply passing a list of lists to groups (without copying a whole anndata object). Will that be sufficient?



sc.set_figure_params(figsize=(20, 20))

cluster_small_multiples(adata=anndata,clust_key='orig.ident',figout=fig_out,palette=sample_palette,size=20)
# sc.pl.umap(tmp, groups=tmp.obs[clust].cat.categories[1:].values, color=adata.obs[clust_key].cat.categories.tolist(), size=size, frameon=frameon, legend_loc=legend_loc, **kwargs)

# sc.pl.umap(anndata,color='orig.ident',palette=sample_palette,save=fig_out,s=50,frameon=False, ncols=4,groups=anndata.obs['orig.ident'].cat.categories.tolist())
fig_out = '_celltypes_umapplots_grouped'+today.strftime("%m%d%y")+'.png'

cluster_small_multiples(adata=anndata,clust_key='Cell_Types_SingleR_scsorter_v2',figout=fig_out,palette=sample_palette,size=20)

sc.tl.leiden(anndata, key_added = "leiden_1.0") # default resolution in 1.0
sc.tl.leiden(anndata, resolution = 0.5, key_added = "leiden_0.5")
sc.tl.leiden(anndata, resolution = 0.7, key_added = "leiden_0.7")
sc.tl.leiden(anndata, resolution = 1.2, key_added = "leiden_1.2")
sc.tl.leiden(anndata, resolution = 1.5, key_added = "leiden_1.5")

results_file = 'GBM.CellState.' + today.strftime("%m%d%y") + 'clustered.h5ad'  # the file that will store the analysis results
anndata.write(results_file)


def make_custom_palette(category_list):
	#rainbow palettes from R in a python list length 50
	rainbow_palette = ["#FF0000","#FF1F00","#FF3D00","#FF5C00","#FF7A00","#FF9900","#FFB800","#FFD600","#FFF500","#EBFF00","#CCFF00","#ADFF00","#8FFF00","#70FF00","#52FF00","#33FF00","#14FF00","#00FF0A","#00FF29","#00FF47","#00FF66","#00FF85","#00FFA3","#00FFC2","#00FFE0","#00FFFF","#00E0FF","#00C2FF","#00A3FF","#0085FF","#0066FF","#0047FF","#0029FF","#000AFF","#1400FF","#3300FF","#5200FF","#7000FF","#8F00FF","#AD00FF","#CC00FF","#EB00FF","#FF00F5","#FF00D6","#FF00B8","#FF0099","#FF007A","#FF005C","#FF003D","#FF001F"]
	thirty_color_pal = ["#fb298a","#2b9b00","#f846ba","#01ca75","#ff84f5","#006a16","#fc3557","#00cb90","#ff5a76","#00c7d9","#f1552b","#008db1","#a33300","#c0aaff","#d6c94a","#324f97","#bc6400","#793882","#9dd48d","#a31925","#018b6a","#ff8670","#185f2a","#8b3358","#5c6100","#ffaa5d","#8f360e","#ffb377","#735800","#9e734d"]
	if(len(category_list) > len(thirty_color_pal)):
		len_sel = len(category_list) - len(thirty_color_pal)
		sel_rainbow = random.sample(rainbow_palette, len_sel)
		sel_pal = thirty_color_pal+sel_rainbow
	else:
		sel_pal = thirty_color_pal
	pal_dict  = {}
	counter = 0
	for i in category_list:
		pal_dict[i] = sel_pal[counter]
		counter = counter+1
	return(pal_dict)


sc.set_figure_params(figsize=(10, 10))
# anndata.obs['Cell_Types_SingleR_scsorter_v2'].cat.categories.tolist()

def plot_clusters_annot(adata,key_word,fig_out,pal,sup_pal):
	tmp = adata.copy()
	leg_title = 'clustering of cells '+str(key_word)
	if(pal==True):
		sc.pl.umap(adata, color=key_word, add_outline=True,legend_loc='on data',
		legend_fontsize=12, legend_fontoutline=2,size=30,legend_fontweight='bold',frameon=False,
		title=leg_title, palette=sup_pal,save=fig_out)
	else:
		sup_pal = make_custom_palette(category_list=adata.obs[key_word].cat.categories.tolist())
		sc.pl.umap(adata, color=key_word, add_outline=True,legend_loc='on data',
		legend_fontsize=12, legend_fontoutline=2,size=30,legend_fontweight='bold',frameon=False,
		title=leg_title,palette=sup_pal,save=fig_out)

for i in ["leiden_0.5","leiden_0.7","leiden_1.0","leiden_1.2","leiden_1.5"]:
	fig_out = '_Clusters_'+str(i)+'_umapplots_annot'+today.strftime("%m%d%y")+'.png'
	plot_clusters_annot(adata=anndata,key_word=i,fig_out=fig_out,pal=False,sup_pal='na')



fig_out = str('_Cell_Types_SingleR_scsorter_v2')+'_umapplots_annot'+today.strftime("%m%d%y")+'.png'

sc.pl.umap(anndata, color=['Cell_Types_SingleR_scsorter_v2'], add_outline=True,legend_loc='right margin',
legend_fontsize=12, legend_fontoutline=2,size=30,legend_fontweight='bold',frameon=False,
title='Cell_Types_SingleR_scsorter_v2', palette=celltype_colorpalette,save=fig_out)


fig_out = str('_orig.ident')+'_umapplots_annot'+today.strftime("%m%d%y")+'.png'

sc.pl.umap(anndata, color=['orig.ident'], add_outline=True,legend_loc='right margin',
legend_fontsize=12, legend_fontoutline=2,size=30,legend_fontweight='bold',frameon=False,
title='orig.ident', palette=sample_palette,save=fig_out)

###Add batch variable to observations

anndata.obs['Batch'] = '5_prime'

old_3prime_batch = list(anndata.obs.index[anndata.obs['orig.ident']=='B152'])+list(anndata.obs.index[anndata.obs['orig.ident']=='B150'])

anndata.obs.loc[old_3prime_batch, 'Batch'] = '3_prime'

batch_palette = {
	'5_prime':"#D95F02",
	'3_prime':"#7570B3"}

fig_out = str('_Batch')+'_umapplots_annot'+today.strftime("%m%d%y")+'.png'

sc.pl.umap(anndata, color=['Batch'], add_outline=True,legend_loc='right margin',
legend_fontsize=12, legend_fontoutline=2,size=30,legend_fontweight='bold',frameon=False,
title='Batch', palette=batch_palette,save=fig_out)

fig_out = str('_Batch_5prime')+'_umapplots_annot'+today.strftime("%m%d%y")+'.png'

sc.pl.umap(anndata, color=['Batch'], add_outline=True,legend_loc='right margin',
legend_fontsize=12, legend_fontoutline=2,size=30,legend_fontweight='bold',frameon=False,
title='Batch',groups=['5_prime'],save=fig_out)


fig_out = str('_Batch_3prime')+'_umapplots_annot'+today.strftime("%m%d%y")+'.png'

sc.pl.umap(anndata, color=['Batch'], add_outline=True,legend_loc='right margin',
legend_fontsize=12, legend_fontoutline=2,size=30,legend_fontweight='bold',frameon=False,
title='Batch',groups=['3_prime'],save=fig_out)



anndata = anndata.read_h5ad('/storage1/fs1/allegra.petti/Active/GBM/CellState/GBM.9/Scanpy_end2end_analysis/GBM.CellState.051722clustered.h5ad')

sc.set_figure_params(scanpy=True, fontsize=20,figsize=(20, 20))

fig_out = str('_n_counts')+'_umapplots_annot'+today.strftime("%m%d%y")+'.png'

sc.pl.umap(anndata, color=['n_counts'], add_outline=True,legend_loc='right margin',
legend_fontsize=12, legend_fontoutline=2,size=10,legend_fontweight='bold',frameon=False,
title='nCount',save=fig_out)

sc.set_figure_params(scanpy=True, fontsize=20,figsize=(30, 30))

for i in ["leiden_0.5","leiden_0.7","leiden_1.0","leiden_1.2","leiden_1.5"]:
	# fig_out = '_Clusters_'+str(i)+'_umapplots_annot'+today.strftime("%m%d%y")+'.png'
	fig_out = str('_n_counts_')+i+'_violinplots.'+today.strftime("%m%d%y")+'.png'
	sc.pl.violin(anndata, ['n_counts'], groupby=i,save=fig_out)
# 	sc.pl.umap(anndata, color=['n_counts'], add_outline=True,legend_loc='right margin',
# legend_fontsize=12, legend_fontoutline=2,size=10,legend_fontweight='bold',frameon=False,
# title='nCount',save=fig_out)

sc.set_figure_params(scanpy=True, fontsize=14,figsize=(25, 25))

for i in ["orig.ident","Batch"]:
	# fig_out = '_Clusters_'+str(i)+'_umapplots_annot'+today.strftime("%m%d%y")+'.png'
	fig_out = str('_n_counts_')+i+'_violinplots.'+today.strftime("%m%d%y")+'.png'
	sc.pl.violin(anndata, ['n_counts'], groupby=i,save=fig_out)

sc.set_figure_params(scanpy=True, fontsize=14,figsize=(20, 20))

for i in anndata.obs['Cell_Types_SingleR_scsorter_v2'].cat.categories.tolist():
	fig_out = str('Cell_Types_SingleR_scsorter_v2_')+i+'_umapplots.'+today.strftime("%m%d%y")+'.png'
	sc.pl.umap(anndata, color=['Cell_Types_SingleR_scsorter_v2'], add_outline=True,legend_loc='right margin',
	legend_fontsize=12, legend_fontoutline=2,size=30,legend_fontweight='bold',frameon=False,
	title=str('Cell type prediction ')+str(i),groups=[i],save=fig_out)

sc.set_figure_params(scanpy=True, fontsize=20,figsize=(20, 20))

#to return to raw anndata
#anndata.raw.to_adata()

for i in anndata.obs['orig.ident'].cat.categories.tolist():
	fig_out = str('orig.ident_')+i+'_umapplots.'+today.strftime("%m%d%y")+'.png'
	sc.pl.umap(anndata, color=['orig.ident'], add_outline=True,legend_loc='right margin',
	legend_fontsize=12, legend_fontoutline=2,size=30,legend_fontweight='bold',frameon=False,
	title=str('orig.ident ')+str(i),groups=[i],save=fig_out)

#####Batch correction approaches#######################

#Trying batch correction methods in scib

out = anndata.copy() #always make copy of anndata object.

# out.raw.to_adata()

sce.pp.harmony_integrate(out, 'Batch')

results_file = 'GBM.CellState.' + today.strftime("%m%d%y") + '.clustered.batch_corrected.harmony.h5ad'  # the file that will store the analysis results
out.write(results_file)

# sc.pl.embedding(out, basis='X_pca_harmony', color=['Batch'],save=fig_out)


adata = out.copy() #always make copy of anndata object.

adata.obsm['X_pca'] = adata.obsm['X_pca_harmony']

sc.pl.umap(adata, color=['Batch'], add_outline=True,legend_loc='right margin',
legend_fontsize=12, legend_fontoutline=2,size=30,legend_fontweight='bold',frameon=False,
title='Batch',save=fig_out)

##### Combat batch correction ######################
ann_combat = anndata.read_h5ad('/storage1/fs1/allegra.petti/Active/GBM/CellState/GBM.9/Scanpy_end2end_analysis/GBM.CellState.batchcorr_combat051922.h5ad')

# out = scib.integration.runScvi(out, 'Batch')

# scvi.data.setup_anndata(out, layer="counts", batch_key = 'Batch')

# seurat_ann  = scib.preprocessing.read_seurat('/storage1/fs1/allegra.petti/Active/GBM/CellState/GBM.9/Redo_PCA_clustering/GBM_cellstate_selgenes.220516.rds')


sc.pp.highly_variable_genes(ann_combat)
print("Highly variable genes: %d"%sum(ann_combat.var.highly_variable))
sc.pl.highly_variable_genes(ann_combat)
sc.pp.pca(ann_combat, n_comps=50, use_highly_variable=True, svd_solver='arpack')
sc.pp.neighbors(ann_combat, n_pcs =50)
sc.tl.umap(ann_combat)

fig_out = str('_Batch_combat_corrected')+'_umapplots_annot'+today.strftime("%m%d%y")+'.png'
sc.pl.umap(ann_combat, color=['Batch'], add_outline=True,legend_loc='right margin',
legend_fontsize=12, legend_fontoutline=2,size=30,legend_fontweight='bold',frameon=False,
title='Batch (Combat batch correction)',save=fig_out)

sc.tl.leiden(ann_combat,resolution = 1.2, key_added = "leiden_1.2_combat")

fig_out = '_Clusters_'+str('leiden_1.2_combat')+'_umapplots_annot'+today.strftime("%m%d%y")+'.png'
plot_clusters_annot(adata=ann_combat,key_word='leiden_1.2_combat',fig_out=fig_out,pal=False,sup_pal='na')


##########################

########## Harmony batch correction ################

adata_harmony = anndata.read_h5ad('/storage1/fs1/allegra.petti/Active/GBM/CellState/GBM.9/Scanpy_end2end_analysis/GBM.CellState.051922.clustered.batch_corrected.harmony.h5ad')
sc.pp.neighbors(adata_harmony, n_pcs =50, use_rep = "X_pca_harmony")
sc.tl.umap(adata_harmony)

fig_out = str('_Batch_harmonycorrected')+'_umapplots_annot'+today.strftime("%m%d%y")+'.png'

# sc.pp.neighbors(out,use_rep='X_pca_harmony',n_neighbors=10, n_pcs=50)
sc.tl.leiden(adata_harmony,resolution = 1.2, key_added = "leiden_1.2_harmony")

fig_out = '_Clusters_'+str('leiden_1.2_harmony')+'_umapplots_annot'+today.strftime("%m%d%y")+'.png'
plot_clusters_annot(adata=adata_harmony,key_word='leiden_1.2_harmony',fig_out=fig_out,pal=False,sup_pal='na')

fig_out = str('_Batch_harmonycorrected')+'_umapplots_annot'+today.strftime("%m%d%y")+'.png'

sc.pl.umap(adata_harmony, color=['Batch'], add_outline=True,legend_loc='right margin',
legend_fontsize=12, legend_fontoutline=2,size=30,legend_fontweight='bold',frameon=False,
title='Batch (harmony batch correction)',save=fig_out)

fig_out = str('orig.ident_harmony_')+'_umapplots.'+today.strftime("%m%d%y")+'.png'
sc.pl.umap(adata_harmony, color=['orig.ident'], add_outline=True,legend_loc='right margin',
legend_fontsize=12, legend_fontoutline=2,size=30,legend_fontweight='bold',frameon=False,
title=str('orig.ident (harmony batch correction)'),save=fig_out)

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors

#Define a nice colour map for gene expression
colors2 = plt.cm.Reds(np.linspace(0, 1, 128))
colors3 = plt.cm.Greys_r(np.linspace(0.7,0.8,20))
colorsComb = np.vstack([colors3, colors2])
mymap = colors.LinearSegmentedColormap.from_list('my_colormap', colorsComb)

fig_out = str('CD45_expression_harmony_')+'_umapplots.'+today.strftime("%m%d%y")+'.png'
sc.pl.umap(adata_harmony, color=['PTPRC'], use_raw=False, color_map=mymap, size = 30,
	title=str('CD45 expression'),save=fig_out)

fig_out = str('Tcellmarkers_harmony_')+'_umapplots.'+today.strftime("%m%d%y")+'.png'
sc.pl.umap(adata_harmony, color=['CD8A','CD8B', 'CD4'], use_raw=False, color_map=mymap, size = 30,
	title=str('T cell markers'),save=fig_out)

fig_out = str('NK/Tmarkers_harmony_')+'_umapplots.'+today.strftime("%m%d%y")+'.png'
sc.pl.umap(adata_harmony, color=['NKG7','CD160', 'GNLY','GZMA'], use_raw=False, color_map=mymap,
	title=str('NK-T/NK markers'),save=fig_out)

fig_out = str('macs_monomarkers_harmony_')+'_umapplots.'+today.strftime("%m%d%y")+'.png'
sc.pl.umap(adata_harmony, color=['FCGR3A','CD14','TREM2'], use_raw=False, color_map=mymap,
	title=str('Monocytes/Macs/Micro-glia'),save=fig_out)

fig_out = str('DC_monocytic_markers_harmony_')+'_umapplots.'+today.strftime("%m%d%y")+'.png'
sc.pl.umap(adata_harmony, color=['CD1C', 'CST3', 'FCER1A'], use_raw=False, color_map=mymap,
	title=str('Monocyte-derived Dendritic cells'),save=fig_out)

fig_out = str('plasmocytoidDC_markers_')+'_umapplots.'+today.strftime("%m%d%y")+'.png'
sc.pl.umap(adata_harmony, color=['GZMB', 'IL3RA'], use_raw=False, color_map=mymap,
title=str('Plasmacytoid dendritic cells'),save=fig_out)

fig_out = str('Bcell_markers_')+'_umapplots.'+today.strftime("%m%d%y")+'.png'

sc.pl.umap(adata_harmony, color=['MS4A1','VPREB3'], use_raw=False, color_map=mymap,
title=str('B cells'),save=fig_out)


##########################
# out = scib.ig.runCombat(out, 'Batch')
