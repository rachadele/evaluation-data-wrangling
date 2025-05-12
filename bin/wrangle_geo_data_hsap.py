#!/bin/python
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
from sklearn.ensemble import RandomForestClassifier
import matplotlib.pyplot as plt
import seaborn as sns
import json
import argparse
import os
import json
from collections import defaultdict
warnings.filterwarnings("ignore")

join_dict = {
    "GSE237718": {
        "meta_key": "SnRNA-seq ID",
        "gemma_key": "sample_name",
        "join_fn": lambda x: str(x).replace("GEX_", "S").replace("R", "")
    }
}



dev_stage_mapping_dict = {
    range(0, 3): "HsapDv_0000083",
    range(2, 6): "HsapDv_0000084",
    range(6, 13): "HsapDv_0000085",
    range(13, 19): "HsapDv_0000086",
    range(19, 44): "HsapDv_0000088",
    range(45, 999): "HsapDv_0000091"
    #range(65, 999): "HsapDv_0000093"
}

study_dev_stage_mapping_dict = {
  "GSE180670":  np.nan 
}

def parse_arguments():
  parser = argparse.ArgumentParser(description="Process .h5ad files.")
  parser.add_argument("--directory", type=str, default="/space/grp/rschwartz/rschwartz/get_gemma_data.nf/study_names_human.txt_author_true_sample_split_true/homo_sapiens", help="Directory containing .h5ad files")
  parser.add_argument("--outdir", type=str, default="/space/grp/rschwartz/rschwartz/evaluation_data_wrangling/pipeline_queries_hsap/sample_subsets", help="Output directory for processed files")
  parser.add_argument("--meta_path", type=str, help="Path to metadata files", default ="/space/grp/rschwartz/rschwartz/evaluation_data_wrangling/meta/GEO_meta")
  if __name__ == "__main__":
    known_args, _ = parser.parse_known_args()
    return known_args
  
  
def read_h5ad_files(directory):
	"""
	Reads all .h5ad files in the specified directory and returns a dictionary of AnnData objects.
	"""
	queries = defaultdict(dict)
	for root, dirs, files in os.walk(directory):
		for file in files:
			if file.endswith(".h5ad"):
				query_name = file.replace(".h5ad", "")
				study_name = file.split("_")[0]
				query = ad.read_h5ad(os.path.join(root, file))
				queries[query_name] = query
	return queries


#def get_sex(query):
	#sex_dict = {
		#"female": "F",
		#"male": "M"
	#}
	#query.obs["sex"] = query.obs["sex"].replace(sex_dict)
	#return query

def get_tissue(query):
  region_dict = {
  "Brain (temporal cortex)": "temporal cortex",
  "adult prefrontal cortex": "prefrontal cortex",
  "brain-resection path (normal tissue) from neurosurgery for epilepsy": "oligodendrocytes in primary culture"}
  query.obs["region"] = query.obs["region"].replace(region_dict)
  return query

def get_specific_meta(meta_path):
  meta_dict= {}
  for root, dirs, files in os.walk(meta_path):
    for file in files:
      if file.endswith(".xlsx"):
       # with open(os.path.join(root, file), "r") as f:
        meta = pd.read_excel(os.path.join(root,file))
        study_name = file.split("_")[0]
        meta_dict[study_name] = meta
  return meta_dict
  
def map_meta(query, meta, join_functions):
  # add value errors to deal with missing keys
  
  gemma_key = join_functions["gemma_key"]
  meta_key = join_functions["meta_key"]
  join_function = join_functions["join_fn"]
  if gemma_key not in query.obs.columns:
    raise ValueError(f"Key {gemma_key} not found in query.obs")
  if meta_key not in meta.columns:
    raise ValueError(f"Key {meta_key} not found in meta")
  query.obs["join_key"] = query.obs[gemma_key].apply(join_function)
  meta["join_key"] = meta[meta_key]
    # this doesn't work because the join key is not in the meta
  new_obs = query.obs.merge(meta, left_on="join_key", right_on="join_key", how="left", suffixes=("", "_y"))
  query.obs = new_obs
  return query


def get_dev_stage(age):
    for age_range, stage in dev_stage_mapping_dict.items():
        if age in age_range:
            return stage
    return "" # For out-of-range values
  
  
def main():

  args = parse_arguments()
  directory = args.directory
  outdir = args.outdir
  os.makedirs(outdir, exist_ok=True)
  if args.meta_path:
    meta_dict = get_specific_meta(meta_path=args.meta_path)
  else:
    meta_dict=None

  # Read all .h5ad files in the specified directory
  queries = read_h5ad_files(directory)
    
  for query_name, query in queries.items():
    study_name = query_name.split("_")[0]
    # if variable is associated with a value
    if meta_dict is not None: 
      meta = meta_dict.get(study_name)
      #if meta_dict.get(study_name) is not None:
      join_functions = join_dict.get(study_name)
      if join_functions is not None and meta is not None:
        query = map_meta(query, meta_dict[study_name], join_functions = join_functions) 
    
    query = get_tissue(query)
    if "Age" in query.obs.columns:
      query.obs["age"] = query.obs["Age"]
      query.obs.drop(columns="Age", inplace=True)
    # if "age" or "Age" in query.obs.columns:
    if "age" in query.obs.columns:
      query.obs["dev_stage"] = query.obs["age"].apply(get_dev_stage)
    else:
      query.obs["dev_stage"] = study_dev_stage_mapping_dict[study_name]

    # drop columns with all NaN values
    query.obs.dropna(axis=1, how='all', inplace=True)
    #query = get_age(query)
    query.write_h5ad(os.path.join(args.outdir, f"{query_name}.h5ad"))
    

  
    
if __name__=="__main__":
	main()


