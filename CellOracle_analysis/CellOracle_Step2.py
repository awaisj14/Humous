 #!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 22 16:44:52 2023

@author: javed
"""

# 0. Import

import os
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

import celloracle as co
co.__version__
# visualization settings
%config InlineBackend.figure_format = 'retina'
%matplotlib inline

plt.rcParams['figure.figsize'] = [6, 4.5]
plt.rcParams["savefig.dpi"] = 300
save_folder = "figures"
os.makedirs(save_folder, exist_ok=True)

adata = sc.read_h5ad("Final_celloracle.h5ad")

df = pd.read_parquet("base_GRN_dataframe_filtered_cicero.parquet")   

base_GRN = df

# Check data
base_GRN.head()

# Instantiate Oracle object
oracle = co.Oracle()


# Check data in anndata
print("Metadata columns :", list(adata.obs.columns))
print("Dimensional reduction: ", list(adata.obsm.keys()))


# In this notebook, we use the unscaled mRNA count for the nput of Oracle object.
adata.X = adata.layers["raw_count"].copy()

# Instantiate Oracle object.
oracle.import_anndata_as_raw_count(adata=adata,
                                   cluster_column_name="celltype",
                                   embedding_name="X_draw_graph_fa")
# You can load TF info dataframe with the following code.
oracle.import_TF_data(TF_info_matrix=base_GRN)

# Alternatively, if you saved the informmation as a dictionary, you can use the code below.
# oracle.import_TF_data(TFdict=TFinfo_dictionary)

#Database was downloaded from: https://cdn.netbiol.org/tflink/download_files/TFLink_Mus_musculus_interactions_All_GMT_proteinName_v1.0.gmt


input_gmt = "/Users/Javed/Documents/Humous/TFLink_Homo_sapiens_interactions_All_GMT_proteinName_v1.0.gmt"
output_csv = "TFLink_Homo_sapiens_interactions_All_GMT_proteinName_v1.0.csv"


# Initialize a list to hold the TF and Target_genes pairs
data = []

# Read and process the GMT file
with open(input_gmt, 'r') as file:
    for line_number, line in enumerate(file, start=1):
        # Strip newline characters and split by tab
        parts = line.strip().split('\t')
        
        # Check if the line has at least three columns (TF, Description, at least one gene)
        if len(parts) < 3:
            print(f"Warning: Line {line_number} in {input_gmt} does not have enough columns. Skipping.")
            continue
        
        tf = parts[0]  # GeneSetName as TF
        # description = parts[1]  # Description (ignored)
        target_genes = parts[2:]  # List of Target_genes
        
        # Append each TF and its Target_genes to the data list
        for gene in target_genes:
            # Optional: Skip empty gene entries
            if gene.strip() == "":
                continue
            data.append({'TF': tf, 'Target_genes': gene})

# Create a DataFrame from the data list
df = pd.DataFrame(data)

# Optional: Remove duplicate entries if any
df.drop_duplicates(inplace=True)

# Save the DataFrame to a CSV file
df.to_csv(output_csv, index=False)

print(f"DataFrame with TF and Target_genes has been saved to {output_csv}")


Paul_15_data = df

# Make dictionary: dictionary key is TF and dictionary value is list of target genes.
TF_to_TG_dictionary = {}

for TF, TGs in zip(Paul_15_data.TF, Paul_15_data.Target_genes):
    # convert target gene to list
    TG_list = TGs.replace(" ", "").split(",")
    # store target gene list in a dictionary
    TF_to_TG_dictionary[TF] = TG_list

# We invert the dictionary above using a utility function in celloracle.
TG_to_TF_dictionary = co.utility.inverse_dictionary(TF_to_TG_dictionary)


# Add TF information
oracle.addTFinfo_dictionary(TG_to_TF_dictionary)
# Perform PCA
oracle.perform_PCA()

# Select important PCs
plt.plot(np.cumsum(oracle.pca.explained_variance_ratio_)[:100])
n_comps = np.where(np.diff(np.diff(np.cumsum(oracle.pca.explained_variance_ratio_))>0.002))[0][0]
plt.axvline(n_comps, c="k")
plt.show()
print(n_comps)
n_comps = min(n_comps, 50)


n_cell = oracle.adata.shape[0]
print(f"cell number is :{n_cell}")
k = int(0.025*n_cell)
print(f"Auto-selected k is :{k}")
oracle.knn_imputation(n_pca_dims=n_comps, k=k, balanced=True, b_sight=k*8,
                      b_maxl=k*4, n_jobs=4)

# Save oracle object.
oracle.to_hdf5("Paul_15_data_v2.celloracle.oracle")
# Load file.
oracle = co.load_hdf5("Paul_15_data_v2.celloracle.oracle")


# Check clustering data
sc.pl.draw_graph(oracle.adata, color="celltype")


# Calculate GRN for each population in "louvain_annot" clustering unit.
# This step may take some time.(~30 minutes)
links = oracle.get_links(cluster_name_for_GRN_unit="celltype", alpha=10,
                         verbose_level=10)
links.links_dict.keys()

# Save Links object.
links.to_hdf5(file_path="links_all.celloracle.links")



## NETWORK PRE PROCESSING
links = co.load_hdf5("links_all.celloracle.links")

links.filter_links(p=0.05, weight="coef_abs", threshold_number=500000)

plt.rcParams["figure.figsize"] = [9, 4.5]
links.plot_degree_distributions(plot_model=True,
                                               #save=f"{save_folder}/degree_distribution/",
                                               )
plt.rcParams["figure.figsize"] = [6, 4.5]
# Calculate network scores.
links.get_network_score()
links.merged_score.head()


# Save Links object.
links.to_hdf5(file_path="links_filtered.celloracle.links")
# You can load files with the following command.
links = co.load_hdf5(file_path="links_filtered.celloracle.links")

# Check cluster name
links.cluster




cluster_name = "RG"
filtered_links_df = links.filtered_links[cluster_name]
filtered_links_df.head()
filtered_links_df.to_csv("GRN_FINAL/RG_GRN_final.csv")

cluster_name = "IPC"
filtered_links_df = links.filtered_links[cluster_name]
filtered_links_df.head()
filtered_links_df.to_csv("GRN_FINAL/IP_GRN_final.csv")

cluster_name = "N"
filtered_links_df = links.filtered_links[cluster_name]
filtered_links_df.head()
filtered_links_df.to_csv("GRN_FINAL/N_GRN_final.csv")











