#!/bin/bash
#SBATCH --job-name=download
#SBATCH -p standard
#SBATCH -A vswarup_lab
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --error=slurm-%J.err
#SBATCH --mem 4G
#SBATCH --array=1-36
#SBATCH --time=4:00:00

# get a file with all of the paths (don't run this in the batch script, just run it before!!)
#cd /dfs7/swaruplab/smorabit/collab/woodlab/cocaine_mouse_2021/
#wget --spider -r --no-parent  https://hts.igb.uci.edu/sudeshnd23041422 2>&1 | grep .fastq | grep -v Removing | tr -s ' ' | cut -d ' ' -f 3 > bin/download_filepaths_April_2023_b1.txt

# output directory
output_dir="/dfs7/swaruplab/smorabit/collab/woodlab/cocaine_mouse_2021/data/April_2023/fastqs/"
cd $output_dir

# download file:
dl_file="/dfs7/swaruplab/smorabit/collab/woodlab/cocaine_mouse_2021/bin/download_filepaths_April_2023_b1.txt"

# set index to the slurm task ID
let index="$SLURM_ARRAY_TASK_ID"

# get the current file to download from the filepaths txt file:
file=$(head -n $index $dl_file | tail -n 1)
echo $file

# check if the file exists:
test_file=$(echo $file | cut -d '/' -f 5)
echo $test_file

if test -f "$test_file"; then
    echo "$test_file exists."
else 
    wget $file
fi
