#!/bin/bash
#SBATCH --job-name=split-by-factors
#SBATCH --partition=Normal
#SBATCH --mem=64G  
#SBATCH --cpus-per-task=4
#SBATCH --nodes=2
# Exit script if any command fails
#set -e

eval "$(conda shell.bash hook)"
conda activate /home/rschwartz/anaconda3/envs/scanpyenv

python /space/grp/rschwartz/rschwartz/evaluation_data_wrangling/bin/split_by_factors.py --query_path /space/grp/rschwartz/rschwartz/evaluation_data_wrangling/pipeline_queries_eric/h5ad/lim.h5ad --query_name lim
python /space/grp/rschwartz/rschwartz/evaluation_data_wrangling/bin/split_by_factors.py --query_path /space/grp/rschwartz/rschwartz/evaluation_data_wrangling/pipeline_queries_eric/h5ad/pineda.h5ad --query_name pineda
python /space/grp/rschwartz/rschwartz/evaluation_data_wrangling/bin/split_by_factors.py --query_path /space/grp/rschwartz/rschwartz/evaluation_data_wrangling/pipeline_queries_eric/h5ad/velmeshev.h5ad --query_name velmeshev 
python /space/grp/rschwartz/rschwartz/evaluation_data_wrangling/bin/split_by_factors.py --query_path /space/grp/rschwartz/rschwartz/evaluation_data_wrangling/pipeline_queries_eric/h5ad/nagy.h5ad --query_name nagy
python /space/grp/rschwartz/rschwartz/evaluation_data_wrangling/bin/split_by_factors.py --query_path /space/grp/rschwartz/rschwartz/evaluation_data_wrangling/pipeline_queries_eric/h5ad/lau.h5ad --query_name lau
python /space/grp/rschwartz/rschwartz/evaluation_data_wrangling/bin/split_by_factors.py --query_path /space/grp/rschwartz/rschwartz/evaluation_data_wrangling/pipeline_queries_eric/h5ad/rosmap.h5ad --query_name rosmap
