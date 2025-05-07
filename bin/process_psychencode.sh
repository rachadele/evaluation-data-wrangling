#!/bin/bash
#SBATCH --job-name=process_nairuz
#SBATCH --partition=Normal
#SBATCH --mem=64G  
#SBATCH --cpus-per-task=4
#SBATCH --nodes=2
# Exit script if any command fails
#set -e

eval "$(conda shell.bash hook)"
conda activate /home/rschwartz/anaconda3/envs/scanpyenv

DATA_DIR=$1

ls "$DATA_DIR" | while read line; do
    full_path="$DATA_DIR/$line"
	echo $full_path
    python /space/grp/rschwartz/rschwartz/evaluation_data_wrangling/bin/wrangle_nairuz_data.py --adata_path $full_path
done
