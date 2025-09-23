#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 16 16:24:35 2023

@author: javed
"""

# Install velocyto from conda using this code if you have M1max MacOS
#conda install -c bioconda velocyto.py

import os
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scvelo as scvelo
import re
import anndata

sc.settings.verbosity = 3             # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.logging.print_versions()
sc.settings.set_figure_params(dpi=600, facecolor='white')

# visualization settings
%config InlineBackend.figure_format = 'retina'
%matplotlib inline

plt.rcParams['figure.figsize'] = [6, 4.5]
plt.rcParams["savefig.dpi"] = 600

adata = sc.read_h5ad("wang_vb_humous_RNA.h5ad") # Preprocessed as the permutations
tre_multi = sc.read_h5ad("updated_data_with_celltype_phase_redone.h5ad") # taken from Bioarchive_version/CellOracle_step1

adata.obs['dataset']
adata.obs['celltype'] = adata.obs['diff_ek']
adata.obs['age_ek']
tre_multi.obs['dataset'] = "TreMulti2021"
tre_multi.obs['age_ek'] = "21"


adata_concat = sc.concat([adata, tre_multi])
adata = adata_concat
sc.pp.filter_genes(adata, min_counts=1)


# Normalize gene expression matrix with total UMI count per cell
sc.pp.normalize_per_cell(adata, key_n_counts='n_counts_all')

# Keep genes that do not contain a dot (.)
canonical_gene_mask = ~adata.var_names.str.contains(r'\.')

# Filter the AnnData object
adata = adata[:, canonical_gene_mask].copy()

# Select top 3000 highly-variable genes
filter_result = sc.pp.filter_genes_dispersion(adata.X,
                                              flavor='cell_ranger',
                                              n_top_genes=3000,
                                              log=False)



# Subset the genes
adata = adata[:, filter_result.gene_subset]

# Renormalize after filtering
sc.pp.normalize_per_cell(adata)


# keep raw cont data before log transformation
adata.raw = adata
adata.layers["raw_count"] = adata.raw.X.copy()

vars = adata.var

adata.obs['dataset']
# Log transformation and scaling
sc.pp.log1p(adata)
sc.pp.scale(adata)

# PCA
sc.tl.pca(adata, svd_solver='arpack')
sc.external.pp.bbknn(adata, batch_key="dataset")


sc.tl.diffmap(adata)
sc.tl.umap(adata)
sc.pl.umap(adata, color="age_ek")
sc.pl.umap(adata, color="dataset")
sc.pl.umap(adata, color="celltype")

adata.obs["celltype"] = adata.obs["celltype"].replace("IP", "IPC")

# PAGA graph construction
sc.tl.paga(adata, groups='celltype')
plt.rcParams["figure.figsize"] = [6, 4.5]
sc.pl.paga(adata)

sc.tl.draw_graph(adata, init_pos='paga', random_state=123)
sc.pl.draw_graph(adata, color='celltype')




G1S_genes_Tirosh = ['MCM5', 'PCNA', 'TYMS', 'FEN1', 'MCM2', 'MCM4', 'RRM1', 'UNG', 'GINS2', 'MCM6', 'CDCA7', 'DTL', 'PRIM1', 'UHRF1', 'MLF1IP', 'HELLS', 'RFC2', 'RPA2', 'NASP', 'RAD51AP1', 'GMNN', 'WDR76', 'SLBP', 'CCNE2', 'UBR7', 'POLD3', 'MSH2', 'ATAD2', 'RAD51', 'RRM2', 'CDC45', 'CDC6', 'EXO1', 'TIPIN', 'DSCC1', 'BLM', 'CASP8AP2', 'USP1', 'CLSPN', 'POLA1', 'CHAF1B', 'BRIP1', 'E2F8']
G2M_genes_Tirosh = ['HMGB2', 'CDK1', 'NUSAP1', 'UBE2C', 'BIRC5', 'TPX2', 'TOP2A', 'NDC80', 'CKS2', 'NUF2', 'CKS1B', 'MKI67', 'TMPO', 'CENPF', 'TACC3', 'FAM64A', 'SMC4', 'CCNB2', 'CKAP2L', 'CKAP2', 'AURKB', 'BUB1', 'KIF11', 'ANP32E', 'TUBB4B', 'GTSE1', 'KIF20B', 'HJURP', 'CDCA3', 'HN1', 'CDC20', 'TTK', 'CDC25C', 'KIF2C', 'RANGAP1', 'NCAPD2', 'DLGAP5', 'CDCA2', 'CDCA8', 'ECT2', 'KIF23', 'HMMR', 'AURKA', 'PSRC1', 'ANLN', 'LBR', 'CKAP5', 'CENPE', 'CTCF', 'NEK2', 'G2E3', 'GAS2L3', 'CBX5', 'CENPA']

cell_cycle_genes =     G1S_genes_Tirosh + G2M_genes_Tirosh
print(len(cell_cycle_genes))

s_genes = G1S_genes_Tirosh
g2m_genes = G2M_genes_Tirosh
cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
print(len(cell_cycle_genes))
print( set(G1S_genes_Tirosh + G2M_genes_Tirosh) - set(cell_cycle_genes) )

print( set(G1S_genes_Tirosh + G2M_genes_Tirosh) &  set(adata.var.index) )

sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
import seaborn as sns

phase_threshold = 0
v1 = adata.obs[ 'S_score']
v2 = adata.obs[ 'G2M_score' ]
str_data_inf = ' Human Multiome '
n_x_subplots = 1
fig = plt.figure(figsize = (20,10) ); c = 0
plt.suptitle(str_data_inf +' '+ str(adata.shape) , fontsize = 20  )

c += 1; fig.add_subplot(1,n_x_subplots ,c)
plt.title('Scanpy phase marks', fontsize = 20 )
v_color =  adata.obs[ 'phase']
ax = sns.scatterplot(x=v1,y=v2, hue = v_color )
# Changing fontsize for the legend: 
plt.setp(ax.get_legend().get_texts(), fontsize=20) # for legend text
plt.setp(ax.get_legend().get_title(), fontsize=20) # for legend title    
plt.xlabel('Scanpy G1S' , fontsize = 20)
plt.ylabel('Scanpy G2M', fontsize = 20 )


plt.show()

if 'RNA_snn_res.0.5' in adata.obs.columns:
    adata.obs.drop(columns=['RNA_snn_res.0.5'], inplace=True)

adata.write_h5ad("All_integrated_CO.h5ad")


import copy
import glob
import time
import os
import shutil
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from tqdm.auto import tqdm



import os, sys, shutil, importlib, glob
from tqdm.notebook import tqdm

#install bedtools2 in your conda environment
#git clone https://github.com/arq5x/bedtools2.git
#cd bedtools2
#make clean; make all
#bin/bedtools --version
#bedtools v2.20.1-4-gb877b35


#import time
import velocyto
import celloracle as co
from celloracle.applications import Pseudotime_calculator
co.__version__

from celloracle import motif_analysis as ma
import celloracle as co
co.__version__

from pybedtools import BedTool



#plt.rcParams["font.family"] = "arial"
plt.rcParams["figure.figsize"] = [5,5]
%config InlineBackend.figure_format = 'retina'
plt.rcParams["savefig.dpi"] = 300

%matplotlib inline




pt = Pseudotime_calculator(adata=adata,
                           obsm_key="X_draw_graph_fa", # Dimensional reduction data name
                           cluster_column_name="celltype" # Clustering data name
                           )


print("Clustering name: ", pt.cluster_column_name)
print("Cluster list", pt.cluster_list)
# Check data
pt.plot_cluster(fontsize=8)



# Here, clusters can be classified into either MEP lineage or GMP lineage

clusters_in_Diff_lineage = ['RG','IPC','N']

# Make a dictionary
lineage_dictionary = {"Lineage_Diff": clusters_in_Diff_lineage}

# Input lineage information into pseudotime object
pt.set_lineage(lineage_dictionary=lineage_dictionary)

# Visualize lineage information
pt.plot_lineages()


import plotly.express as px
import plotly.io as pio
pio.renderers.default = "browser"
try:
    import plotly.express as px
    def plot(adata, embedding_key, cluster_column_name):
        embedding = adata.obsm[embedding_key]
        df = pd.DataFrame(embedding, columns=["x", "y"])
        df["cluster"] = adata.obs[cluster_column_name].values
        df["label"] = adata.obs.index.values
        fig = px.scatter(df, x="x", y="y", hover_name=df["label"], color="cluster")
        fig.show()

    plot(adata=pt.adata,
         embedding_key=pt.obsm_key,
         cluster_column_name=pt.cluster_column_name)
except:
    print("Plotly not found in your environment. Did you install plotly? Please read the instruction above.")

# Estimated root cell name for each lineage
root_cells = {"Lineage_Diff": "hF_WK8.5cer.ctx_2019_ACGTACATCTCCACTG"}
pt.set_root_cells(root_cells=root_cells)


pt.plot_root_cells()

"umap" in pt.adata.obsm

pt.get_pseudotime_per_each_lineage()

pt.plot_pseudotime(cmap="rainbow")




adata.obs = pt.adata.obs

pt = adata.obs["Pseudotime"].copy()

# Clip extreme values above 99th percentile
upper = pt.quantile(0.99)
pt_clipped = pt.clip(upper=upper)

# Rescale only within clipped range
pt_min = pt_clipped.min()
pt_max = pt_clipped.max()
adata.obs["Pseudotime_rescaled"] = (pt_clipped - pt_min) / (pt_max - pt_min)

# Check
adata.obs["Pseudotime_rescaled"].describe()

from sklearn.preprocessing import QuantileTransformer

qt = QuantileTransformer(output_distribution="uniform", random_state=0)
adata.obs["Pseudotime_quantile"] = qt.fit_transform(pt.values.reshape(-1, 1))

sc.pl.umap(adata, color="Pseudotime_rescaled", cmap="viridis")
sc.pl.umap(adata, color="Pseudotime_quantile", cmap="viridis")

adata.obs["Pseudotime"] =adata.obs["Pseudotime_quantile"]

sc.pl.draw_graph(adata, color='Pseudotime_quantile')



adata.write_h5ad("Final_celloracle.h5ad")

adata = sc.read_h5ad("Final_celloracle.h5ad")


# Load scATAC-seq peak list.
peaks = pd.read_csv("wang_all_peaks.csv", index_col=0)
peaks.replace('-','_', inplace=True,regex=True)
peaks = peaks.x.values

# Load Cicero coaccessibility scores.
cicero_connections = pd.read_csv("wang_cicero_connections.csv", index_col=0)
cicero_connections.head()









ref_genome = "hg38"

genome_installation = ma.is_genome_installed(ref_genome=ref_genome)
print(ref_genome, "installation: ", genome_installation)
if not genome_installation:
    import genomepy
    genomepy.install_genome(name=ref_genome, provider="UCSC")
else:
    print(ref_genome, "is installed.")
tss_annotated = ma.get_tss_info(peak_str_list=peaks, ref_genome="hg38")

# Check results
tss_annotated.tail()

integrated = ma.integrate_tss_peak_with_cicero(tss_peak=tss_annotated,
                                               cicero_connections=cicero_connections)
print(integrated.shape)
integrated.head()

peak = integrated[integrated.coaccess >= 0.8]
peak = peak[["peak_id", "gene_short_name"]].reset_index(drop=True)
print(peak.shape)
peak.head()

peak.to_csv("processed_peak_file_wang.csv")
# Load both dataframes
peaks_tre = pd.read_csv("processed_peak_file.csv", index_col=0)
peak = pd.read_csv("processed_peak_file_wang.csv", index_col=0)

# Combine them vertically
combined_peaks = pd.concat([peaks_tre, peak], ignore_index=True)

# Optional: drop duplicates based on peak_id
combined_peaks = combined_peaks.drop_duplicates(subset="peak_id")

# Save to CSV
combined_peaks.to_csv("combined_peaks.csv", index=False)





# Load annotated peak data.
peaks = pd.read_csv("combined_peaks.csv")
peaks.head()

def decompose_chrstr(peak_str):
    """
    Args:
        peak_str (str): peak_str. e.g. 'chr1_3094484_3095479'

    Returns:
        tuple: chromosome name, start position, end position
    """

    *chr_, start, end = peak_str.split("_")
    chr_ = "_".join(chr_)
    return chr_, start, end

from genomepy import Genome

def check_peak_format(peaks_df, ref_genome):
    """
    Check peak format.
     (1) Check chromosome name.
     (2) Check peak size (length) and remove sort DNA sequences (<5bp)

    """

    df = peaks_df.copy()

    n_peaks_before = df.shape[0]

    # Decompose peaks and make df
    decomposed = [decompose_chrstr(peak_str) for peak_str in df["peak_id"]]
    df_decomposed = pd.DataFrame(np.array(decomposed), index=peaks_df.index)
    df_decomposed.columns = ["chr", "start", "end"]
    df_decomposed["start"] = df_decomposed["start"].astype(int)
    df_decomposed["end"] = df_decomposed["end"].astype(int)

    # Load genome data
    genome_data = Genome(ref_genome)
    all_chr_list = list(genome_data.keys())


    # DNA length check
    lengths = np.abs(df_decomposed["end"] - df_decomposed["start"])


    # Filter peaks with invalid chromosome name
    n_threshold = 5
    df = df[(lengths >= n_threshold) & df_decomposed.chr.isin(all_chr_list)]

    # DNA length check
    lengths = np.abs(df_decomposed["end"] - df_decomposed["start"])

    # Data counting
    n_invalid_length = len(lengths[lengths < n_threshold])
    n_peaks_invalid_chr = n_peaks_before - df_decomposed.chr.isin(all_chr_list).sum()
    n_peaks_after = df.shape[0]


    #
    print("Peaks before filtering: ", n_peaks_before)
    print("Peaks with invalid chr_name: ", n_peaks_invalid_chr)
    print("Peaks with invalid length: ", n_invalid_length)
    print("Peaks after filtering: ", n_peaks_after)

    return df


peaks = check_peak_format(peaks, ref_genome)

# Instantiate TFinfo object
tfi = ma.TFinfo(peak_data_frame=peaks,
                ref_genome=ref_genome)

# Scan motifs. !!CAUTION!! This step may take several hours if you have many peaks!
tfi.scan(fpr=0.02,
         motifs=None,  # If you enter None, default motifs will be loaded.
         verbose=True)

# Save tfinfo object
tfi.to_hdf5(file_path="test1.celloracle.tfinfo")
# Check motif scan results
tfi.scanned_df.head()


# Reset filtering
tfi.reset_filtering()

# Do filtering
tfi.filter_motifs_by_score(threshold=1)

# Format post-filtering results.
tfi.make_TFinfo_dataframe_and_dictionary(verbose=True)


df = tfi.to_dataframe()
df.head()


# Save result as a dataframe
df = tfi.to_dataframe()
df.to_parquet("base_GRN_dataframe_filtered_cicero.parquet")




