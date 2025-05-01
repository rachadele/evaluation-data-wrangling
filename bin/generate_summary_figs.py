#!/user/bin/python3

from pathlib import Path
import random
import os
import sys
import scanpy as sc
import numpy as np
import pandas as pd
import anndata as ad
import cellxgene_census
import scvi
from scipy.sparse import csr_matrix
import warnings
import cellxgene_census
import cellxgene_census.experimental
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
import json
import argparse
import os
import json
import yaml

import sys
sys.path.append("/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/bin")
import adata_functions
from adata_functions import *
from adata_functions import relabel

# Function to parse command line arguments
def parse_arguments():
    parser = argparse.ArgumentParser(description="Plot cell type level statistics based on organism, census version, and label mapping file.")
    parser.add_argument('--organism', type=str, default='mus_musculus', help='Organism name (e.g., homo_sapiens)')
    parser.add_argument('--census_version', type=str, default='2024-07-01', help='Census version (e.g., 2024-07-01)')
    parser.add_argument('--relabel_path', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/census_map_mouse.tsv")
    parser.add_argument('--ref_collections', type=str, nargs = '+', default = [
       "A taxonomy of transcriptomic cell types across the isocortex and hippocampal formation",
        "An integrated transcriptomic and epigenomic atlas of mouse primary motor cortex cell types",
        "Adult mouse cortical cell taxonomy revealed by single cell transcriptomics"
    ]) 
    parser.add_argument('--split_column', type=str, default="tissue")
    parser.add_argument('--ref_keys', type=str, nargs="+", default=["subclass","class","family","global"])
    parser.add_argument('--source_data_dir', type=str, default="/space/grp/rschwartz/rschwartz/get_gemma_data.nf/forebrain_only/mus_musculus")
    parser.add_argument('--seed', type=int, default=42)
    parser.add_argument('--relabel_dir', type=str, default="/space/grp/rschwartz/rschwartz/nextflow_eval_pipeline/meta/relabel_mus_musculus")
    if __name__ == "__main__":
        known_args, _ = parser.parse_known_args()
        return known_args

def plot_celltype_contributions(combined_obs, outpath, columns_to_plot=['subclass', 'class', 'family'], split="dataset_title"):
    """
    Creates a bar plot showing the contribution of each query to various categories 
    (subclass, class, family) and saves the figure.

    Parameters:
    - combined_obs: DataFrame containing the combined observations from all queries.
    - outdir: Directory path where the figure will be saved.
    - columns_to_plot: List of columns (e.g., ['subclass', 'class', 'family']) to plot contributions for.
    """
    # Set up subplots (one for each column in columns_to_plot)
    fig, axes = plt.subplots(1, len(columns_to_plot), figsize=(len(columns_to_plot) * 10, 10), sharey=True)
    
    # Loop through the columns and create a plot for each
    for i, col in enumerate(columns_to_plot):
        # Count occurrences of query per category
        data_counts = combined_obs.groupby([col, split]).size().unstack(fill_value=0)
        
        # Drop categories with all zero counts
        # this isn't working for "split" column
        
        data_counts = data_counts.loc[:, (data_counts != 0).any(axis=0)] 
        #drop columns with all zero counts
        data_counts = data_counts.loc[:, (data_counts != 0).any(axis=0)] 
        # Normalize to proportions
        data_counts = data_counts.div(data_counts.sum(axis=1), axis=0)
        
        # Plot
        ax = axes[i]
        data_counts.plot(kind="bar", stacked=True, ax=ax, colormap="tab10")
        ax.tick_params(axis='x', labelsize=22, rotation=90)
        ax.set_xlabel("")
        ax.set_ylabel("Proportion")
        ax.set_title(f"{col.capitalize()}")
        if i == len(columns_to_plot) - 1:
            ax.legend(title=split, bbox_to_anchor=(1, 1))
        else:
            ax.legend().remove()
        ax.tick_params(axis='x', rotation=90)

    # Save figure
    plt.savefig(os.path.join(outpath), dpi=300, bbox_inches="tight")


def main():
     # Parse command line arguments
    args = parse_arguments()

    plt.rcParams['font.size'] = 35         # Set default font size
    # Set organism and census_version from arguments
    organism = args.organism
    census_version = args.census_version
    relabel_path = args.relabel_path
    ref_collections=args.ref_collections
    split_column=args.split_column
    ref_keys=args.ref_keys
    SEED = args.seed
    source_data_dir = args.source_data_dir
    outdir=f"/space/grp/rschwartz/rschwartz/evaluation_data_wrangling/census_summary/{organism}"
    os.makedirs(outdir, exist_ok=True)
    census = cellxgene_census.open_soma(census_version=census_version)
    dataset_info = census.get("census_info").get("datasets").read().concat().to_pandas()
    brain_obs = cellxgene_census.get_obs(census, organism,
        value_filter=(
            "tissue_general == 'brain' and "
            "is_primary_data == True and "
            "disease == 'normal' "
        ))

    brain_obs = brain_obs.merge(dataset_info, on="dataset_id", suffixes=(None,"_y"))
    brain_obs.drop(columns=['soma_joinid_y'], inplace=True)
    brain_obs_filtered = brain_obs[brain_obs['collection_name'].isin(ref_collections)] 
    celltype_breakdown=brain_obs_filtered[["cell_type","dataset_title","collection_name"]].value_counts().reset_index()
    celltype_breakdown.to_csv(os.path.join(outdir, "celltype_breakdown.tsv"),sep="\t", index=False)

    
    relabel_df = pd.read_csv(relabel_path, sep='\t')  # Adjust the separator as needed
    # Take the first column as the join key
    join_key = relabel_df.columns[0]
    # Ensure the join_key is in both the AnnData object and the relabel DataFrame
    if join_key not in brain_obs_filtered.columns:
        raise ValueError(f"{join_key} not found in AnnData object observations.")
    if join_key not in relabel_df.columns:
        raise ValueError(f"{join_key} not found in relabel DataFrame.")
    # Perform the left join to update the metadata
    brain_obs_filtered = brain_obs_filtered.merge(relabel_df, on=join_key, how='left', suffixes=(None, "_y"))
  
    # Plot
    plot_celltype_contributions(brain_obs_filtered, outpath=os.path.join(outdir, "dataset_ref_celltype_contribution.png"), columns_to_plot=ref_keys, split="dataset_title")
    plot_celltype_contributions(brain_obs_filtered, outpath=os.path.join(outdir, "tissue_ref_celltype_contribution.png"), columns_to_plot=ref_keys, split="tissue")
    queries= {}
    for filepath in os.listdir(source_data_dir):
        if filepath.endswith(".h5ad"):
            name = os.path.splitext(filepath)[0]
            filepath = os.path.join(source_data_dir, filepath)
            queries[name] = ad.read_h5ad(filepath)
            
    relabel_dir = args.relabel_dir
    
    new_queries = {}
    for name in queries.keys():
        name_new = name.split("_")[0]
        relabel_path = os.path.join(relabel_dir, f"{name_new}_relabel.tsv")
        if os.path.exists(relabel_path):
            new_queries[name] = relabel(queries[name], relabel_path)
            new_queries[name].obs["dataset_title"] = name_new
            
    common_columns = set.intersection(*[set(adata.obs.columns) for adata in new_queries.values()])

    # Combine metadata from all queries using only common columns
    combined_obs = pd.concat(
        [adata.obs[list(common_columns)].assign(query=name) for name, adata in new_queries.items()],
        ignore_index=True
    )
    
    # make sure all columns are strings
    combined_obs = combined_obs.astype(str)
    combined_obs = combined_obs.fillna("None")

   
# Example usage:
    plot_celltype_contributions(combined_obs, outpath=os.path.join(outdir, "query_study_cell_type_contribution.png"), columns_to_plot=ref_keys, split="dataset_title")
    plot_celltype_contributions(combined_obs, outpath=os.path.join(outdir, "query_region_cell_type_contribution.png"), columns_to_plot=ref_keys, split="region") 
    plot_celltype_contributions(combined_obs, outpath=os.path.join(outdir, "query_region_cell_type_contribution.png"), columns_to_plot=ref_keys, split="region") 
    if organism == "homo_sapiens":
        plot_celltype_contributions(combined_obs, outpath=os.path.join(outdir, "query_region_cell_type_contribution.png"), columns_to_plot=ref_keys, split="disease") 
        plot_celltype_contributions(combined_obs, outpath=os.path.join(outdir, "query_region_cell_type_contribution.png"), columns_to_plot=ref_keys, split="dev_stage") 
 
    #if organism == "mus_musculus":
     #   plot_celltype_contributions(combined_obs, outpath=os.path.join(outdir, "query_region_cell_type_contribution.png"), columns_to_plot=ref_keys, split="genotype")
      #  plot_celltype_contributions(combined_obs, outpath=os.path.join(outdir, "query_region_cell_type_contribution.png"), columns_to_plot=ref_keys, split="treatment")
       # plot_celltype_contributions(combined_obs, outpath=os.path.join(outdir, "query_region_cell_type_contribution.png"), columns_to_plot=ref_keys, split="strain")

        
        
if __name__ == "__main__":  
    main()
