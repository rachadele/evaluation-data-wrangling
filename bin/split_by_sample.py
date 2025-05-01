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
  parser.add_argument('--query_path', type=str, default="/space/grp/rschwartz/rschwartz/evaluation_data_wrangling/pipeline_queries_eric/lim.h5ad")
  parser.add_argument('--outdir', default="/space/grp/rschwartz/rschwartz/evaluation_data_wrangling/pipeline_queries_eric/sample_subsets")
  
  if __name__ == "__main__":
    known_args, _ = parser.parse_known_args()
    return known_args
  
def subset_by_sample(adata, query_name, sample_key, outdir):
    subsets = {}
    for sample in adata.obs[sample_key].unique():
        subsets[sample] = adata[adata.obs[sample_key] == sample]
        subsets[sample].write_h5ad(os.path.join(outdir,f"{query_name}_{sample}.h5ad"))

def main():
  args = parse_arguments()
  query_path = args.query_path
  query_name = query_path.split("/")[-1].replace(".h5ad", "")  # More robust file name extraction
  adata = ad.read_h5ad(query_path)
  outdir = args.outdir

  sample_keys = {
  "lau": "sample",
  "pineda": "Sample_ID",
  "lim": "case_num",
  "velmeshev": "sample",
  "rosmap": "individualID", 
  "nagy": "batch" # old = in barcode, extracted, no sample info
  }
  
  os.makedirs(outdir, exist_ok=True)
  
  
  sample_key = sample_keys[query_name]
  subset_by_sample(adata, query_name, sample_key, outdir)
    
       

if __name__ == "__main__":
  main()