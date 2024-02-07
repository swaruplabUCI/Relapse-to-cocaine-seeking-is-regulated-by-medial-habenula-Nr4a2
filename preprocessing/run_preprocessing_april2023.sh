#!/bin/bash
#SBATCH --job-name=process
#SBATCH -p standard
#SBATCH -A vswarup_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --error=slurm-%J.err
#SBATCH --mem 16G
#SBATCH --array=0-7
#SBATCH --time=1:00:00

source ~/.bashrc
conda activate scrublet

data_dir="/dfs7/swaruplab/smorabit/collab/woodlab/cocaine_mouse_2021/data/April_2023/cellranger/"

samples=($(ls $data_dir))
let index="$SLURM_ARRAY_TASK_ID"
sample=${samples[$index]}

# create output directory
mkdir $data_dir$sample"/outs/scanpy_preprocessing"

# process for cellranger data:
python /dfs7/swaruplab/smorabit/collab/woodlab/cocaine_mouse_2021/bin/scanpy_preprocessing.py \
  --infile $data_dir$sample"/outs/filtered_feature_bc_matrix.h5" \
  --outdir $data_dir$sample"/outs/scanpy_preprocessing/" \
  --mt 5 \
  --maxcount 20000 \
  --mincount 500 \
  --leidenres 0.5 \
  --npcs 30 \
  --neighbors 20 \
  --mindist 0.25

# run if we have run cellbender
if test -f $data_dir$sample/outs/cellbender_filtered.h5; then

  mkdir $data_dir$sample"/outs/scanpy_preprocessing_cellbender"

  python /dfs7/swaruplab/smorabit/collab/woodlab/cocaine_mouse_2021/bin/scanpy_preprocessing.py \
    --infile $data_dir$sample"/outs/cellbender_filtered.h5" \
    --outdir $data_dir$sample"/outs/scanpy_preprocessing_cellbender/" \
    --mt 5 \
    --maxcount 20000 \
    --mincount 500 \
    --leidenres 0.5 \
    --npcs 30 \
    --neighbors 20 \
    --mindist 0.25
fi
