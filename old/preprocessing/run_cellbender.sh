#!/bin/bash

# activate appropriate env
source ~/.bashrc
conda activate CellBender

cd /dfs3b/swaruplab/smorabit/collab/woodlab/cocaine_mouse_2021/test_cellbender/

# get a list of all samples:
cellranger_dir="/dfs3b/swaruplab/smorabit/collab/woodlab/cocaine_mouse_2021/data/cellranger/"

# loop through each sample
for sample in $cellranger_dir/*; do
  echo $sample

  # set input and output file
  infile=$sample/outs/raw_feature_bc_matrix.h5
  outfile=$sample/outs/cellbender.h5

  # run cellbender
  cellbender remove-background \
       --input $infile \
       --output $outfile \
       --expected-cells 10000 \
       --total-droplets-included 25000 \
       --epochs 150 \
       --cuda

done;
