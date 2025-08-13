import os
import pandas as pd
import numpy as np
from collections import defaultdict
import seaborn as sns
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.patches as mpatches

# Paths and setup
parent_dir = "/home/h/hm435/Desktop/Clone/Clonesig_Results2"
sns.set(context="talk")

# Containers: per-sample signature counts split by clonal/subclonal
clonal_counts = defaultdict(lambda: defaultdict(int))
subclonal_counts = defaultdict(lambda: defaultdict(int))

# Parse all subfolders and aggregate counts
for subdir in os.listdir(parent_dir):
    sub_path = os.path.join(parent_dir, subdir)
    if not os.path.isdir(sub_path):
        continue

    for file in os.listdir(sub_path):
        if not file.endswith(".csv"):
            continue

        df = pd.read_csv(os.path.join(sub_path, file))
        if df.empty or "sample_id" not in df.columns:
            continue

        sample_id = df["sample_id"].iloc[0]

        for _, row in df.iterrows():
            sig = f"SBS{int(row['signature']) + 1}"  # signatures are 1-indexed in labels
            clone_id = int(row["clonesig_clone_reassignment"])
            if clone_id == 0:
                clonal_counts[sample_id][sig] += 1
            else:
                subclonal_counts[sample_id][sig] += 1

# Collect full signature list (sorted numerically) and sample list
all_signatures = sorted(
    {sig for d in [clonal_counts, subclonal_counts] for s in d.values() for sig in s},
    key=lambda x: int(x.replace("SBS", ""))
)
samples = sorted(set(clonal_counts.keys()) | set(subclonal_counts.keys()))

# convert nested dict -> DataFrame [samples x signatures]
def build_df(counts_dict):
    data = [[counts_dict[sample].get(sig, 0) for sig in all_signatures] for sample in samples]
    return pd.DataFrame(data, index=samples, columns=all_signatures)

clonal_df = build_df(clonal_counts)
subclonal_df = build_df(subclonal_counts)

# Row-wise normalize (per sample); 0 if a row is all zeros
clonal_df = clonal_df.div(clonal_df.max(axis=1), axis=0).fillna(0)
subclonal_df = subclonal_df.div(subclonal_df.max(axis=1), axis=0).fillna(0)

# Combine clonal (left block) and subclonal (right block)
combined_df = pd.concat([clonal_df, subclonal_df], axis=1)
combined_df.columns = ["clonal"] * len(clonal_df.columns) + ["subclonal"] * len(subclonal_df.columns)
signature_labels = all_signatures * 2  # labels for the two blocks

# Clustermap heatmap
g = sns.clustermap(
    combined_df,
    cmap="Greens",
    metric="euclidean",
    figsize=(14, 12),
    xticklabels=False,
    yticklabels=False,
    col_cluster=False,
    cbar_pos=None
)

# Annotate axes and add split line between blocks
ax = g.ax_heatmap
num_sigs = len(all_signatures)
ax.axvline(x=num_sigs, color="black", linewidth=2)
ax.set_xticks(np.arange(2 * num_sigs) + 0.5)
ax.set_xticklabels(signature_labels, rotation=90, fontsize=22)
ax.text(num_sigs / 2, -15, "Clonal", ha="center", va="center", fontsize=23)
ax.text(num_sigs + num_sigs / 2, -15, "Subclonal", ha="center", va="center", fontsize=23)
ax.set_xlabel("COSMIC SBS Signatures", fontsize=22, labelpad=20)
ax.set_ylabel("Samples", fontsize=14, labelpad=22)

# Manual colorbar above the row dendrogram
ax_dendro = g.ax_row_dendrogram
ax_dendro.text(
    0.055, 0.5, "Samples", fontsize=18, ha="center", va="center",
    rotation=90, transform=ax_dendro.transAxes
)
cax = inset_axes(
    ax_dendro, width="100%", height="5%", loc="upper center",
    bbox_to_anchor=(0.05, 0.12, 1, 1),
    bbox_transform=ax_dendro.transAxes, borderpad=0
)
norm = plt.Normalize(vmin=0, vmax=1)
sm = plt.cm.ScalarMappable(cmap="Greens", norm=norm)
sm.set_array([])
cbar = plt.colorbar(sm, cax=cax, orientation="horizontal")
cbar.ax.text(1.05, 0.48, "Intensity", transform=cbar.ax.transAxes,
             va="center", ha="left", fontsize=20)

# Title and layout
g.fig.suptitle(
    "Clonal vs Subclonal Signatures Extracted with CloneSig (N=277)",
    fontsize=24, x=0.57, y=0.93
)
g.fig.subplots_adjust(left=0.01, bottom=0.2, right=0.99, top=0.958)

# Save heatmap
g.savefig("clonal_subclonal_heatmap_277.png", dpi=300, bbox_inches="tight")
plt.show()

# Subclonal minus clonal mean intensity (per signature)
signature_deltas = subclonal_df.mean() - clonal_df.mean()
signature_deltas = signature_deltas.sort_values(ascending=False)
print("Subclonal - Clonal Mean Signature Intensity Differences:")
print(signature_deltas)

# Figure 2: bar plot of deltas
fig, ax = plt.subplots(figsize=(10, 6))
colors = ["green" if v > 0 else "red" for v in signature_deltas]
signature_deltas.plot(kind="bar", color=colors, ax=ax)

ax.tick_params(axis="x", labelsize=12)
ax.axhline(0, color="black", linewidth=0.8)
ax.set_ylabel("Subclonal - Clonal Mean Intensity", fontsize=15)
ax.set_title("Signature Activity Change (Subclonal - Clonal) Extracted with CloneSig", fontsize=16, y=1.05)

green_patch = mpatches.Patch(color="green", label=">0 higher in subclonal")
red_patch = mpatches.Patch(color="red", label="<0 higher in clonal")
ax.legend(handles=[green_patch, red_patch],
          loc="upper center", bbox_to_anchor=(0.5, -0.25),
          ncol=2, fontsize=12, frameon=False)

fig.subplots_adjust(left=0.1, bottom=0.28, top=0.88, right=0.74)

# Optional: annotate N and (if available) significant p-values
ax.text(0.98, 0.95, "N = 277", transform=ax.transAxes, fontsize=13,
        va="top", ha="right",
        bbox=dict(boxstyle="square,pad=0.3", edgecolor="black", facecolor="white"))

fig.savefig("clonal_subclonal_change_277.png", dpi=300)
plt.show()
