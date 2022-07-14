import scanpy as sc
import scrublet as scr
import numpy as np
from matplotlib import pyplot as plt
import argparse
import os
sc.settings.set_figure_params(dpi=500, dpi_save=1000, figsize=(5,5), facecolor='white')

# set up arguments
parser = argparse.ArgumentParser()
parser.add_argument("--infile", help="path to the input file", type=str)
parser.add_argument("--outdir", help="directory to store outputs", type=str)
parser.add_argument("--mt", help="Filter out cells with greater than this percent of mitochondrial reads", type=float, default=5.0)
parser.add_argument("--maxcount", help="Filter out cells with greater than or equal to this number of UMI", type=float, default=None)
parser.add_argument("--mincount", help="Filter out cells with fewer than or equal to this number of UMI", type=float, default=0)
parser.add_argument("--maxgene", help="Filter out cells with greater than or equal to this number of genes detected", type=float, default=None)
parser.add_argument("--mingene", help="Filter out cells with fewer than or equal to this number of genes detected", type=float, default=0)
parser.add_argument("--leidenres", help="Resolution parameter for Leiden clustering", type=float, default=1.0)
parser.add_argument("--npcs", help="The number of Principal Components to include for UMAP", type=int, default=30)
parser.add_argument("--neighbors", help="The number of nearest neighbors to include for UMAP", type=int, default=20)
parser.add_argument("--mindist", help="Minimum distane between two points in UMAP space", type=float, default=0.2)

# parse arguments
args = parser.parse_args()

infile = args.infile
outdir = args.outdir
mt_thresh = args.mt
max_count = args.maxcount
min_count = args.mincount
max_gene = args.maxgene
min_gene = args.mingene
leiden_res = args.leidenres
n_pcs = args.npcs
n_neighbors = args.neighbors
min_dist = args.mindist

################################################################################
# sanitize inputs
################################################################################

if outdir[-1] != '/':
    outdir += '/'

################################################################################
# create output directories:
################################################################################

fig_dir = '{}figures/'.format(outdir)
data_dir = '{}data/'.format(outdir)

for d in fig_dir, data_dir:
    if not os.path.exists(d):
        os.makedirs(d)

sc.settings.figdir = fig_dir

################################################################################
# Load data and predict doublets
################################################################################

adata = sc.read_10x_h5(infile, genome='mm10')
adata.var_names_make_unique()

scrub = scr.Scrublet(adata.X)
adata.obs['doublet_scores'], predicted_doublets = scrub.scrub_doublets()
adata.obs['doublets'] = np.where(predicted_doublets== False, 'Singlet', 'Doublet')

# plot the scrublet histogram
ax = scrub.plot_histogram()
plt.savefig('{}scrublet_histogram.pdf'.format(fig_dir))

# save the resulting adata object:
adata.write('{}adata_unfiltered.h5ad'.format(data_dir))


################################################################################
# scanpy processing
################################################################################

# percent mito
adata.var['mt'] = adata.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
sc.pl.violin(adata,['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], inner='box', size=0,  groupby='doublets', multi_panel=True, save='_qc.pdf')

# write the metadata as a csv:
adata.obs.to_csv('{}cell_meta_unprocessed.-=csv'.format(data_dir))

# cell filtering:
adata = adata[adata.obs.pct_counts_mt <= mt_thresh]

# filter by max UMI:
if max_count is not None:
    adata = adata[adata.obs.total_counts <= max_count]

# filter by min UMI:
adata = adata[adata.obs.total_counts >= min_count]

# filter by max gene:
if max_gene is not None:
    adata = adata[adata.obs.n_genes_by_counts <= max_gene]

# filter by min gene:
adata = adata[adata.obs.n_genes_by_counts >= min_gene]

# plot the QC violin plot after filtering:
sc.pl.violin(adata,['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], inner='box', size=0,  groupby='doublets', multi_panel=True, save='_qc_post_filter.pdf')

# pre-processing
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, svd_solver='arpack')

################################################################################
# Non-linear dim reduction + clustering
################################################################################

sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs, metric='cosine')

sc.tl.umap(adata, min_dist=min_dist)
sc.tl.leiden(adata, resolution=leiden_res)

# compute paga clusters
sc.tl.paga(adata)
sc.pl.paga(
    adata, threshold=0.2, title='', edge_width_scale=0.5,
    frameon=False, save=True
)

# re-compute UMAP using paga as initial position
sc.tl.umap(adata, min_dist=min_dist, init_pos='paga')

################################################################################
# plotting
################################################################################

# plot clusters
sc.pl.umap(adata, color=['leiden'],frameon=False, legend_loc='on data', legend_fontoutline=1, legend_fontsize=9, add_outline=True, title='', save='_leiden.pdf')

# plot doublets
sc.pl.umap(adata, color=['doublet_scores', 'doublets', 'pct_counts_mt', 'n_genes_by_counts', 'total_counts'],frameon=False, ncols=3, save='_qc.pdf')

#compute dendrogram for the dotplot:
sc.tl.dendrogram(adata, 'leiden')

# marker genes from dotplot:
marker_dict = {
    'MF-ODC' : ['Plp1', 'Sox10'],
    'Mature ODC' : ['Mog', 'Mobp'],
    'Neuron' : ['Stmn2', 'Thy1', 'Nr4a2', 'Snap25'],
    'MHb_Neuron': ['Tac2', 'Slc17a7'],
    'LHb_Neuron': ['Pcdh10'],
    'Microglia': ['Tmem199', 'Csf1r'],
    'Ependyma':['Foxj1', 'Fam216b'],
    'Astrocyte' : ['Ntsr2', 'Slc4a4'],
    'OPC' : ['Pdgfra'],
    'Mural' : ['Cldn5'],
    'Endothelial': ['Flt1'],
    'Pericyte': ['Pdgfrb'],
    'VLMC': ['Vtn']
}

sc.pl.dotplot(adata, marker_dict, groupby='leiden', figsize=(12,6), standard_scale='var', dendrogram=True, save='markers.pdf')
sc.pl.matrixplot(adata, marker_dict, groupby='leiden', standard_scale='var', dendrogram=True, save='markers.pdf')
adata.write('{}adata_processed.h5ad'.format(data_dir))
