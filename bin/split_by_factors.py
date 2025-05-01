#!/user/bin/python3

from pathlib import Path
import os
import sys
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import warnings
import re
from pathlib import Path
import itertools
import argparse

# Function to parse command line arguments
def parse_arguments():
  parser = argparse.ArgumentParser(description="Download model file based on organism, census version, and tree file.")
  parser.add_argument('--query_path', type=str, default="/space/grp/rschwartz/rschwartz/evaluation_data_wrangling/pipeline_queries_eric/h5ad/velmeshev.h5ad")
  parser.add_argument('--factors', type=str, nargs="+", default=["sex", "disease", "dev_stage", "region"])
  parser.add_argument('--query_name', type=str, default="velmeshev")
  parser.add_argument('--outdir', type=str, default="/space/grp/rschwartz/rschwartz/evaluation_data_wrangling/pipeline_queries_all_files/factor_subsets")
  
  if __name__ == "__main__":
    known_args, _ = parser.parse_known_args()
    return known_args
  

def main():
  args = parse_arguments()
  query = ad.read_h5ad(args.query_path)
  query_name = args.query_name
  factors = args.factors
  outdir = args.outdir
  os.makedirs(outdir, exist_ok=True)
  # i need to split on a grid for each factor so they're not confounded
  # Check if each factor exists in the query
  factors = [factor for factor in factors if factor in query.obs.columns]
  if factors == []:
    raise ValueError("None of the factors provided are present in the query.")
  # Generate all combinations of factor levels
  factor_combinations = itertools.product(*(query.obs[factor].unique() for factor in factors))
  
  # Iterate over all combinations
  for combination in factor_combinations:
    print(combination)
    # Start with the full query object
    query_subset = query
    # Apply filtering condition progressively
    for i, factor in enumerate(factors):
      print(i)
      print(factor)
      # Filter query by the current factor and its corresponding value in the combination
      query_subset = query_subset[query_subset.obs[factor] == combination[i]]

      # Check if the query_subset is empty before proceeding
      if query_subset.n_obs == 0:
          break
      # Only write the subset if it's not empty
    if query_subset.n_obs > 0:
      combination_str = "_".join(map(str, combination)).replace(" ", "_")
      print(combination_str)
    
      query_subset.write_h5ad(os.path.join(outdir,f"{query_name}_{combination_str}.h5ad"))

if __name__ == "__main__":
  main()