#!/bin/bash

# activate appropriate env
source ~/.bashrc
conda activate cellbender

cd /dfs7/swaruplab/smorabit/collab/woodlab/cocaine_mouse_2021/data/April_2023/

# get a list of all samples:
cellranger_dir="/dfs7/swaruplab/smorabit/collab/woodlab/cocaine_mouse_2021/data/April_2023/cellranger/"

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
