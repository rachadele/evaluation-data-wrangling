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



queries = {}
for root,dirs, files in os.walk("/space/grp/rschwartz/rschwartz/evaluation_data_wrangling/pipeline_queries_eric/h5ad"):
    for file in files:
        if file.endswith(".h5ad"):
            query_name = file.replace(".h5ad","")
            query=ad.read_h5ad(os.path.join(root,file))
            query
            queries[query_name] = query
os.chdir("/space/grp/rschwartz/rschwartz/evaluation_data_wrangling/pipeline_queries_eric/h5ad")

nagy = queries["nagy"]
nagy.obs["meta"] = nagy.obs.index
nagy.obs["disease"] = nagy.obs["meta"].str.split("_").apply(lambda x: next((i for i in x if i in ["Suicide", "Control"]), None))
nagy.obs["disease"] = nagy.obs["disease"].str.replace("Suicide","MDD")
nagy.obs["batch"] = nagy.obs["meta"].str.split("_").apply(lambda x: next((i for i in x if re.match(r"^B[1-9]", i)), None))
nagy.obs.drop(columns="meta", inplace=True)
nagy.obs["dev_stage"] = "HsapDv_0000091"
nagy.obs["region"] = "dorsolateral prefrontal cortex"
nagy.obs["sex"] = "M"
nagy.obs["age"] = 38.71
nagy.write_h5ad("nagy.h5ad")


velmeshev = queries["velmeshev"]
velmeshev.obs["region"] = velmeshev.obs["region"].str.replace("PFC","prefrontal cortex")
velmeshev.obs["region"] = velmeshev.obs["region"].str.replace("ACC","anterior cingulate cortex")
velmeshev.obs["dev_stage"] = velmeshev.obs["age"].apply(get_dev_stage)
velmeshev.obs.rename(columns={"diagnosis":"disease"})
velmeshev.write_h5ad("velmeshev.h5ad")

lim = queries["lim"]
lim_meta = pd.read_excel("/space/grp/rschwartz/rschwartz/evaluation_data_wrangling/pipeline_queries_eric/meta/lim_41467_2022_35388_MOESM8_ESM.xlsx")
lim_meta["id"] = lim_meta["Donor"].str.replace("T-","")
lim_meta["id"] = lim_meta["id"] + "_" + lim_meta["Region"]
lim.obs["id"] = lim.obs["case_num"].str.replace(r'[a-zA-Z]', '', regex=True)
lim.obs["id"] = lim.obs["id"].astype(str) + "_" + lim.obs["Location"].astype(str)
lim.obs = lim.obs.merge(lim_meta, how="left", on="id")
lim.obs.rename(columns={"Gender":"sex"}, inplace=True)
lim.obs["disease"] = lim.obs["Condition"].str.replace("Con","Control")
lim.obs.rename(columns={"Age":"age"}, inplace=True)

# missing info for case_num = C5382Cd
# age = 62
# sex = m
# disease = Control

lim.obs["age"] = np.where(lim.obs["case_num"] == "C5382Cd", 62, lim.obs["age"])
lim.obs["sex"] = np.where(lim.obs["case_num"] == "C5382Cd", "M", lim.obs["sex"])
lim.obs["disease"] = np.where(lim.obs["case_num"] == "C5382Cd", "Control", lim.obs["disease"])
                          
lim.obs["dev_stage"] = lim.obs["age"].apply(get_dev_stage)
lim.obs["region"] = lim.obs["Location"].str.replace("Cingulate","cingulate cortex")
lim.obs["region"] = lim.obs["region"].str.replace("Caudate","caudate nucleus")
lim.obs["region"] = lim.obs["region"].str.replace("Accumbens","nucleus accumbens")
lim.write_h5ad("lim.h5ad")


lau = queries["lau"]
lau_meta = pd.read_excel("/space/grp/rschwartz/rschwartz/evaluation_data_wrangling/pipeline_queries_eric/meta/lau_pnas.2008762117.sd01.xlsx")
lau_meta.rename(columns={"ID":"active.ident", "SEX":"sex"}, inplace=True)
lau.obs=lau.obs.merge(lau_meta, left_on="active.ident", right_on="active.ident")
lau.obs["disease"] = lau.obs["active.ident"].str[:2].replace("NC","Control")
lau.obs.drop(columns="AGE-RELATED PLAQUE SCORE", inplace=True)
lau.obs.rename(columns={"AGE":"age"}, inplace=True)
lau.obs["dev_stage"] = lau.obs["age"].apply(get_dev_stage)
lau.obs["region"] = "prefrontal cortex"
lau.write_h5ad("lau.h5ad")

pineda = queries["pineda"]
pineda_meta = pd.read_excel("/space/grp/rschwartz/rschwartz/evaluation_data_wrangling/pipeline_queries_eric/meta/pineda_mmc1.xlsx")
pineda_meta = pineda_meta[ :-4]
# make pineda_meta donor an int
pineda_meta["Donor"] = pineda_meta["Donor"].astype(int)
# make pineda mayo ID a string
pineda.obs["Mayo_ID"] = pineda.obs["Mayo_ID"].astype(int)
# donors don't match between pineda and pineda_meta
#pineda.obs = pineda.obs.merge(pineda_meta, left_on="Mayo_ID", right_on="Donor", suffixes=("","_y"))
# drop columns ending in y
pineda.obs = pineda.obs.loc[:,~pineda.obs.columns.str.endswith("_y")]
pineda.obs.rename(columns={"Sex":"sex"}, inplace=True)
pineda.obs["sex"]=pineda.obs["sex"].str.replace("FALSE","F")
pineda.obs["disease"] = pineda.obs["Condition"].str.replace("PN","Control")
pineda.obs["age"] = pineda.obs["Age of death (Y)"]
pineda.obs["dev_stage"] = "HsapDv_0000091" 
pineda.obs["region"] = "primary motor cortex"
pineda.write_h5ad("pineda.h5ad")

rosmap = queries["rosmap"]
rosmap_meta = {}
for root, dir, file in os.walk("/space/scratch/ericchu/r_cache/041_CH4_FINAL/data/rosmap"):
    for f in file:
        if f.endswith(".csv"):
            rosmap_meta[f] = pd.read_csv(os.path.join(root,f))
        if f.endswith("xlsx"):  
            rosmap_meta[f] = pd.read_excel(os.path.join(root,f), header=1)
        if f.endswith("tsv"):
            rosmap_meta[f] = pd.read_csv(os.path.join(root,f), sep="\t")
            
good_meta = rosmap_meta["ROSMAP_Clinical_2019-05_v3.csv"]
rosmap.obs = rosmap.obs.merge(good_meta, left_on="projid", right_on="projid", suffixes=("","_y"))

columns_to_drop = [column for column in rosmap.obs.columns if column.endswith("_y")]
rosmap.obs.drop(columns=columns_to_drop, inplace=True)
rosmap.obs["batch"] = "batch_" + rosmap.obs["projid"].astype(str)
rosmap.obs["sex"] = rosmap.obs["msex"].replace({1:"M", 0:"F"})
rosmap.obs["age"] = pd.to_numeric(rosmap.obs["age_death"], errors="coerce")
rosmap.obs["age"] = rosmap.obs["age"].fillna(90).astype(int)
rosmap.obs["disease"] = rosmap.obs["cogdx"].replace({1:"Control", 2:"Control", 3:"Control", 4:"AD", 5:"AD", 6:"DEM" }).astype(str)
rosmap.obs["dev_stage"] = rosmap.obs["age"].apply(get_dev_stage)
rosmap.obs["region"] = "prefrontal cortex"
rosmap.write_h5ad("rosmap.h5ad")
