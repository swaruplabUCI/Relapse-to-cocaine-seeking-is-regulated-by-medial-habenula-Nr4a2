
# script to run TF network analysis in parallel 
library(optparse)

option_list = list(
  make_option(
    c('--seurat'), type='character', default=NULL,
    help='Seurat .rds file', metavar='character'
  ),
  make_option(
    c('--wgcnaname'), type='character', default=NULL,
    help='Name of the hdWGCNA experiment in the Seurat object, which already has metacells calculated.', metavar='character'
  ),
  make_option(
    c('--outdir'), type='character', default='./',
    help='Directory to place output files', metavar='character'
  ),
  make_option(
    c('--groupby'), type='character', default=NULL,
    help='Column in the seurat_obj@meta.data slot indicating how to group cells (clusters, celltypes, etc).', metavar='character'
  ),
  make_option(
    c('--index'), type='numeric', default=NULL,
    help='SLURM task array number goes here, selects which group to process on this task'
  ),
  make_option(
    c('--name'), type='character', default='DEGs',
    help='Name to append to output files.'
  ),
  # make_option(
  #   c('--multigroupby'), type='character', default=NULL,
  #   help='Column in the seurat_obj@meta.data slot indicating how to subgroup cells (diagnosis, condition, sample, etc).', metavar='character'
  # ),
  # make_option(
  #     c('--multigroupname'), type='character', default=NULL,
  #     help='comma separated list of groups that are in seurat_obj@meta.data[[multigroupby]] to help us subset the data before running the TF network.', metavar='character'
  # ),
  make_option(
      c('--objective'), type='character', default='reg:squarederror',
      help='xgboost parameter option for objective.', metavar='character'
  ),
  make_option(
      c('--maxdepth'), type='numeric', default=1,
      help='xgboost parameter option for tree depth.', metavar='character'
  ),
  make_option(
      c('--eta'), type='numeric', default=0.1,
      help='xgboost parameter option for eta.', metavar='character'
  ),
  make_option(
      c('--nthread'), type='numeric', default=16,
      help='xgboost parameter option for number of parallel threads.', metavar='character'
  ),
  make_option(
      c('--alpha'), type='numeric', default=0.1,
      help='xgboost parameter option for alpha.', metavar='character'
  ),
  make_option(
      c('--nestimators'), type='numeric', default=50,
      help='xgboost parameter option for the number of estimators.', metavar='character'
  ),
  make_option(
      c('--nfold'), type='numeric', default=5,
      help='xgboost parameter option for the number of folds in cross validation.', metavar='character'
  )
)

# parse arguments
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)

print(opt)

library(Seurat);
library(tidyverse);
library(xgboost);
library(tictoc);

# co-expression network analysis packages:
library(WGCNA);
library(hdWGCNA);
library(igraph);

# set random seed for reproducibility
set.seed(12345)

# load the TF functions (later just include them in hdWGCNA):
source('/dfs7/swaruplab/smorabit/collab/woodlab/cocaine_mouse_2021/Nurr2c_vs_GFP/revision/bin/TF_functions.R')

# load Seurat object:
seurat_obj <- readRDS(opt$seurat)

print('Data loaded successfully!')

#---------------------------------------------------------------------------
# unpack options
#---------------------------------------------------------------------------

wgcna_name <- opt$wgcnaname
out_dir <- opt$outdir
group.by <- opt$groupby
# group_name <- opt$group_name 
# multi.group.by <- opt$multigroupby 
# multi_group_name <- opt$multigroupname

# xgboost options
model_params <- list(
  objective = opt$objective,
  max_depth = opt$maxdepth,
  eta = opt$eta,
  nthread=opt$nthread,
  alpha=opt$alpha,
  n_estimators=opt$nestimators
)
nfold <- opt$nfold

#---------------------------------------------------------------------------
# setup data for the TF network analysis
#---------------------------------------------------------------------------

print('Setting up expression matrix...')

# get a list of all groups:
cell_groups <- seurat_obj@meta.data[,opt$groupby] %>% unique %>% as.character

# get current cluster based on the index:
group_name <- cell_groups[opt$index]

# get motif data
motif_matrix <- GetMotifMatrix(seurat_obj)
motif_df <- GetMotifs(seurat_obj)
motif_df <- subset(motif_df, gene_name %in% rownames(seurat_obj))

# select genes that are in the motif matrix & in the seurat object
wgcna_genes <- GetWGCNAGenes(seurat_obj, wgcna_name=wgcna_name)
genes_use <- wgcna_genes[wgcna_genes %in% rownames(motif_matrix)]
genes_use <- unique(c(genes_use, motif_df$gene_name))
seurat_obj <- SetWGCNAGenes(seurat_obj, genes_use, wgcna_name=wgcna_name)

# running an updated version of SetDatExpr that's not pushed to the repo yet
seurat_obj <- SetDatExpr(
  seurat_obj,
  group.by = group.by,
  group_name = group_name,
  # multi.group.by = multi.group.by,
  # multi_group_name = multi_group_name,
  wgcna_name=wgcna_name
)


#---------------------------------------------------------------------------
# Run the TF network analysis
#---------------------------------------------------------------------------
print('Constructing TF Network...')

net_results <- ConstructTFNetwork(
  seurat_obj,
  model_params=model_params,
  nfold=nfold,
  wgcna_name=wgcna_name
)
importance_df <- net_results[[1]]
eval_df <- net_results[[2]]

# save the results:
write.csv(importance_df, file=paste0(out_dir, '/', opt$name, '_', group_name, '_importance.csv'), quote=FALSE, row.names=FALSE)
write.csv(eval_df, file=paste0(out_dir, '/', opt$name, '_', group_name, '_eval.csv'), quote=FALSE, row.names=FALSE)
