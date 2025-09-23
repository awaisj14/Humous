#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Simulate and visualize CellOracle shift vectors in memory-efficient batches.
Each batch simulates on a subset of the full dataset.
Only RG cells are extracted and visualized (quiver + pie chart), averaged across batches.

NEW:
- Optional random downsampling AFTER binning to limit plotted arrows and avoid overcrowding.
- Safe embedding retrieval; minor robustness tweaks.
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
import scanpy as sc
import celloracle as co
import gc  # Garbage collection

# === CONFIG ===
goi = "JUNB"
perturbation_type = "KO"  # or "KO"
n_batches = 5            # number of batches
arrow_scale = 1250
bin_size = 5.0            # size for spatial binning before averaging
save_folder = f"{goi}_{perturbation_type}_batched"
os.makedirs(save_folder, exist_ok=True)

# Reproducibility
batch_random_seed = 123           # affects which cells go into which batch
plot_random_seed = 42             # affects which binned arrows get plotted
max_plot_arrows = 5000            # None to disable; otherwise plot at most this many binned arrows

# === Matplotlib and PDF settings ===
mpl.rcParams['pdf.fonttype'] = 42
mpl.rcParams['ps.fonttype'] = 42
plt.rcParams["figure.figsize"] = [6, 6]
plt.rcParams["savefig.dpi"] = 600

# === Direction color mapping ===
direction_colors = {
    "G1":  "#dd1c77",  # North
    "G2M": "#c994c7",  # South
    "IPC": "#86c9c2",  # West
    "RG":  "#e2c683",  # East
}

def get_custom_color(angle):
    """Map vector angle (radians) to a direction color."""
    if np.pi/4 <= angle < 3*np.pi/4:
        return direction_colors["G1"]     # North
    elif -3*np.pi/4 <= angle < -np.pi/4:
        return direction_colors["G2M"]    # South
    elif angle >= 3*np.pi/4 or angle < -3*np.pi/4:
        return direction_colors["IPC"]    # West
    else:
        return direction_colors["RG"]     # East

def hex_to_rgb(hex_color):
    return mcolors.to_rgba(hex_color)[:3]

# === Load full oracle (for counting) and links ===
oracle_base = co.load_hdf5("Paul_15_data_v2.celloracle.oracle")
links = co.load_hdf5("links_filtered.celloracle.links")
links.filter_links()

# === Dynamically determine batch sizes ===
total_ncells = oracle_base.adata.n_obs
cells_per_batch = int(np.ceil(total_ncells / n_batches))
print(f"🧠 Total cells: {total_ncells}")
print(f"📦 Dividing into {n_batches} batches of ~{cells_per_batch} cells each.")

rng_batches = np.random.default_rng(batch_random_seed)
all_indices = rng_batches.permutation(total_ncells)
batches = np.array_split(all_indices, n_batches)

for i, b in enumerate(batches):
    print(f"  - Batch {i+1}: {len(b)} cells")

# Free base oracle memory
del oracle_base
gc.collect()

# === Store RG shifts and coords across batches ===
rg_coords_all = []
rg_shifts_all = []

for i, batch_idx in enumerate(batches):
    print(f"\n🚀 Running batch {i+1} of {n_batches} (n = {len(batch_idx)} cells)")

    # Load and subset oracle
    oracle = co.load_hdf5("Paul_15_data_v2.celloracle.oracle")
    # Subset by positional indices (assumes consistent ordering across loads)
    oracle.adata = oracle.adata[batch_idx].copy()

    # Build TF dictionary and simulate
    oracle.get_cluster_specific_TFdict_from_Links(links)
    oracle.fit_GRN_for_simulation(alpha=10, use_cluster_specific_TFdict=True)

    # Get perturbation value
    if perturbation_type == "KO":
        perturb_value = 0.0
    elif perturbation_type == "OE":
        expr = sc.get.obs_df(oracle.adata, keys=[goi], layer="imputed_count")[goi]
        perturb_value = float(expr.max())
    else:
        raise ValueError("perturbation_type must be 'KO' or 'OE'")

    # Run simulation
    oracle.simulate_shift(perturb_condition={goi: perturb_value}, n_propagation=3)
    oracle.estimate_transition_prob(n_neighbors=200, knn_random=True, sampled_fraction=1)
    oracle.calculate_embedding_shift(sigma_corr=0.05)

    # Extract embedding coordinates safely
    if "X_draw_graph_fa" in oracle.adata.obsm:
        coords = oracle.adata.obsm["X_draw_graph_fa"]
    else:
        coords = oracle.adata.obsm[oracle.embedding_name]

    shifts = oracle.delta_embedding

    # Extract RG cells
    if "celltype" not in oracle.adata.obs.columns:
        raise KeyError("`celltype` column not found in adata.obs.")
    rg_mask = (oracle.adata.obs["celltype"] == "RG").to_numpy()

    rg_coords_all.append(coords[rg_mask])
    rg_shifts_all.append(shifts[rg_mask])

    # Clean up
    del oracle
    gc.collect()

# === Combine all RG coordinates and shifts ===
if len(rg_coords_all) == 0 or sum(arr.shape[0] for arr in rg_coords_all) == 0:
    raise RuntimeError("No RG cells found across batches; nothing to plot.")

rg_coords = np.vstack(rg_coords_all)
rg_shifts = np.vstack(rg_shifts_all)
print(f"\n📊 Combined RG vectors: {rg_coords.shape[0]} cells")

# === Bin and average vectors ===
x_bins = np.round(rg_coords[:, 0] / bin_size)
y_bins = np.round(rg_coords[:, 1] / bin_size)
binned_keys = list(zip(x_bins, y_bins))

bin_dict = {}
for key, coord, shift in zip(binned_keys, rg_coords, rg_shifts):
    if key not in bin_dict:
        bin_dict[key] = {"coords": [], "shifts": []}
    bin_dict[key]["coords"].append(coord)
    bin_dict[key]["shifts"].append(shift)

avg_coords = np.array([np.mean(v["coords"], axis=0) for v in bin_dict.values()])
avg_shifts = np.array([np.mean(v["shifts"], axis=0) for v in bin_dict.values()])

print(f"🧮 Bins: {len(bin_dict)}  → averaged vectors: {avg_coords.shape[0]}")

# === Optional random downsampling for plotting (to avoid overcrowding) ===
num_bins = avg_coords.shape[0]
if (max_plot_arrows is not None) and (num_bins > max_plot_arrows):
    rng_plot = np.random.default_rng(plot_random_seed)
    keep_idx = rng_plot.choice(num_bins, size=max_plot_arrows, replace=False)
    avg_coords_plot = avg_coords[keep_idx]
    avg_shifts_plot = avg_shifts[keep_idx]
    print(f"✂️ Downsampling plotted arrows: {num_bins} → {len(keep_idx)}")
else:
    avg_coords_plot = avg_coords
    avg_shifts_plot = avg_shifts
    print(f"✂️ Downsampling skipped (plotting {num_bins} arrows).")

# === Quiver plot (using downsampled arrays) ===
U, V = avg_shifts_plot[:, 0], avg_shifts_plot[:, 1]
angles = np.arctan2(V, U)
arrow_colors = [get_custom_color(a) for a in angles]

fig, ax = plt.subplots(figsize=(6, 6))
ax.quiver(avg_coords_plot[:, 0], avg_coords_plot[:, 1], U * arrow_scale, V * arrow_scale,
          color=arrow_colors, pivot="tail", angles="xy", scale_units="xy", scale=1)
ax.set_title(f"Simulated RG shifts (batched): {goi} {perturbation_type}")
ax.set_xlabel("FA embedding 1")
ax.set_ylabel("FA embedding 2")
ax.set_aspect("equal", "box")
plt.savefig(os.path.join(save_folder, f"batched_quiver_RG_only_{goi}_{perturbation_type}.pdf"))
plt.close(fig)

# === Pie chart of shift directions (of plotted subset) ===
direction_labels = {
    direction_colors["G1"]: "G1 (North)",
    direction_colors["G2M"]: "G2M (South)",
    direction_colors["IPC"]: "IPC (West)",
    direction_colors["RG"]: "RG (East)"
}
color_counts = {v: 0 for v in direction_colors.values()}
for c in arrow_colors:
    color_counts[c] += 1

labels = [direction_labels[c] for c in color_counts.keys()]
sizes = list(color_counts.values())
pie_colors = list(color_counts.keys())

fig_pie, ax_pie = plt.subplots(figsize=(5, 5))
ax_pie.pie(sizes, labels=labels, colors=pie_colors, autopct="%1.1f%%",
           startangle=140, textprops={"fontsize": 8})
ax_pie.set_title("RG Shift Directions (batched simulation; plotted subset)", fontsize=10)
plt.savefig(os.path.join(save_folder, f"batched_pie_RG_only_{goi}_{perturbation_type}.pdf"))
plt.close(fig_pie)

print(f"\n✅ Saved figures to: {save_folder}")
