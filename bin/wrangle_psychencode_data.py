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

dev_stage_mapping_dict = {
    range(0, 3): "HsapDv_0000083",
    range(2, 6): "HsapDv_0000084",
    range(6, 13): "HsapDv_0000085",
    range(13, 19): "HsapDv_0000086",
    range(19, 44): "HsapDv_0000088",
    range(45, 999): "HsapDv_0000091",
    #range(65, 999): "HsapDv_0000093"
}

def get_dev_stage(age):
    for age_range, stage in dev_stage_mapping_dict.items():
        if age in age_range:
            return stage
    return ""  # For out-of-range values


# Function to parse command line arguments
def parse_arguments():
  parser = argparse.ArgumentParser(description="Download model file based on organism, census version, and tree file.")
  parser.add_argument('--adata_path', type=str, default="/space/grp/rschwartz/rschwartz/get_gemma_data.nf/nairuz_data/homo_sapiens/Velmeshev_et_al_1203427.h5ad")
  parser.add_argument('--outdir', default="/space/grp/rschwartz/rschwartz/evaluation_data_wrangling/pipeline_queries_hsap/sample_subsets")
  
  if __name__ == "__main__":
    known_args, _ = parser.parse_known_args()
    return known_args

def main():
    args = parse_arguments()
    adata_path = args.adata_path
    adata = sc.read_h5ad(adata_path)
    outdir =args.outdir
    
    filename = adata_path.split("/")[-1]
    print(filename)

    adata.obs.rename(columns={"Age_death": "age", "Biological Sex": "sex", "Disorder":"disease"}, inplace=True)
    adata.obs["region"] = "prefrontal cortex"
    adata.obs["dev_stage"] = adata.obs["age"].astype(int).apply(get_dev_stage)
    
    adata.write_h5ad(os.path.join(outdir, filename))
    
if __name__ == "__main__":
    main()