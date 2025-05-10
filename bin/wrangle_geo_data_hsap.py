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
    "GSE237718": "NA",
}


dev_stage_mapping_dict = {
    range(0, 3): "HsapDv_0000083",
    range(2, 6): "HsapDv_0000084",
    range(6, 13): "HsapDv_0000085",
    range(13, 19): "HsapDv_0000086",
    range(19, 44): "HsapDv_0000088",
    range(45, 999): "HsapDv_0000091",
    #range(65, 999): "HsapDv_0000093"
}

study_dev_stage_mapping_dict = {
  "GSE237718": "HsapDv_0000091",
  "GSE180670": "unknown"
}



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
  "brain-resection path (normal region) from neurosurgery for epilepsy": "oligodencrocytes in primary culture"}
  query.obs["region"] = query.obs["region"].replace(region_dict)
  return query

def get_specific_meta(meta_path):
  meta_dict= {}
  for root, dirs, files in os.walk(meta_path):
    for file in files:
      if file.endswith(".xlsx"):
        with open(os.path.join(root, file), "r") as f:
          meta = pd.read_excel(f)
        query_name = file.split("_")[0]
        meta_dict[query_name] = meta
  return meta
  
def map_meta(query, meta, join_key):
  new_obs = query.obs.merge(meta, left_on=join_key, right_on=join_key, how="left", suffixes=("", "_y"))
  query.obs = new_obs
  return query


def get_dev_stage(age):
    for age_range, stage in dev_stage_mapping_dict.items():
        if age in age_range:
            return stage
    return "" # For out-of-range values
  
  
def main():
  parser = argparse.ArgumentParser(description="Process .h5ad files.")
  parser.add_argument("--directory", type=str, default="/space/grp/rschwartz/rschwartz/get_gemma_data.nf/hsap_mmus_forebrain_small_samples/homo_sapiens", help="Directory containing .h5ad files")
  parser.add_argument("--outdir", type=str, default="/space/grp/rschwartz/rschwartz/evaluation_data_wrangling/pipeline_queries_hsap/test", help="Output directory for processed files")
  parser.add_argument("--meta_path", type=str, help="Path to metadata files")
  args = parser.parse_args()
  directory = args.directory
  outdir = args.outdir
  os.makedirs(outdir, exist_ok=True)
  if args.meta_path:
    meta_dict = get_specific_meta(args.meta_path)
  else:
    meta_dict=None
  # Read all .h5ad files in the specified directory
  queries = read_h5ad_files(directory)
    
  for query_name, query in queries.items():
    study_name = query_name.split("_")[0]
    # if variable is associated with a value
    if meta_dict is not None and meta_dict.get(query_name) is None: 
      query = map_meta(query, meta_dict[query_name], join_key = join_dict[study_name])
      
    query = get_tissue(query)
    #query = get_sex(query)
    if "Age" in query.obs.columns:
      query.obs["age"] = query.obs["Age"].astype(float)
    if "age" in query.obs.columns:
      query.obs["dev_stage"] = query.obs["age"].apply(get_dev_stage)
    else:
      query.obs["dev_stage"] = study_dev_stage_mapping_dict[study_name]

    #query = get_age(query)
    query.write_h5ad(os.path.join(args.outdir, f"{query_name}.h5ad"))
    

  
    
if __name__=="__main__":
	main()


