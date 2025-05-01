#!/bin/bash
#SBATCH --job-name=split-by-sample
#SBATCH --partition=Normal
#SBATCH --mem=64G  
#SBATCH --cpus-per-task=4
#SBATCH --nodes=2
# Exit script if any command fails
#set -e

eval "$(conda shell.bash hook)"
conda activate /home/rschwartz/anaconda3/envs/scanpyenv

python /space/grp/rschwartz/rschwartz/evaluation_data_wrangling/bin/split_by_sample.py --query_path /space/grp/rschwartz/rschwartz/evaluation_data_wrangling/pipeline_queries_eric/h5ad/lim.h5ad
python /space/grp/rschwartz/rschwartz/evaluation_data_wrangling/bin/split_by_sample.py --query_path /space/grp/rschwartz/rschwartz/evaluation_data_wrangling/pipeline_queries_eric/h5ad/pineda.h5ad
python /space/grp/rschwartz/rschwartz/evaluation_data_wrangling/bin/split_by_sample.py --query_path /space/grp/rschwartz/rschwartz/evaluation_data_wrangling/pipeline_queries_eric/h5ad/velmeshev.h5ad
python /space/grp/rschwartz/rschwartz/evaluation_data_wrangling/bin/split_by_sample.py --query_path /space/grp/rschwartz/rschwartz/evaluation_data_wrangling/pipeline_queries_eric/h5ad/nagy.h5ad
python /space/grp/rschwartz/rschwartz/evaluation_data_wrangling/bin/split_by_sample.py --query_path /space/grp/rschwartz/rschwartz/evaluation_data_wrangling/pipeline_queries_eric/h5ad/lau.h5ad
python /space/grp/rschwartz/rschwartz/evaluation_data_wrangling/bin/split_by_sample.py --query_path /space/grp/rschwartz/rschwartz/evaluation_data_wrangling/pipeline_queries_eric/h5ad/rosmap.h5ad
