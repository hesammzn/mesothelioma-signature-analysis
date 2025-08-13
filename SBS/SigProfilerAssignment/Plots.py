# ============================== Imports ==============================
import os
import itertools
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.gridspec as gridspec
from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.stats import wilcoxon, spearmanr


# ==================== Figure 1: Clonal vs Subclonal heatmaps ====================
# Load activities (tab-sep; first column is sample index)
df_clo = pd.read_csv(
    "C:/Users/hesam/Desktop/irp5/277/Sig/SigAssign_3.4_2/Clonal/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",
    sep="\t", index_col=0
)
df_sub = pd.read_csv(
    "C:/Users/hesam/Desktop/irp5/277/Sig/SigAssign_3.4_2/Subclonal/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",
    sep="\t", index_col=0
)

# Keep columns (signatures) with total counts > 350
df_sub = df_sub.loc[:, df_sub.sum(axis=0) > 350]
df_clo = df_clo.loc[:, df_clo.sum(axis=0) > 350]

# Row-wise max normalization (per sample) → [0,1]
df_sub_norm = df_sub.div(df_sub.max(axis=1), axis=0).fillna(0)
df_clo_norm = df_clo.div(df_clo.max(axis=1), axis=0).fillna(0)

# Ward clustering to reorder rows
Z_clo = linkage(df_clo_norm, method='ward')
Z_sub = linkage(df_sub_norm, method='ward')
ordered_clo = df_clo_norm.iloc[leaves_list(Z_clo)]
ordered_sub = df_sub_norm.iloc[leaves_list(Z_sub)]

# Layout: two heatmaps + shared colorbar column
fig = plt.figure(figsize=(14, 6))
gs = gridspec.GridSpec(1, 3, width_ratios=[5, 5, 0.3], wspace=0.03)

# Left panel: Clonal heatmap (no colorbar)
ax_clo = plt.subplot(gs[0])
sns.heatmap(ordered_clo, cmap="Blues", ax=ax_clo,
            xticklabels=True, yticklabels=False, vmin=0, vmax=1, cbar=False)
ax_clo.set_title("Clonal Signatures", fontsize=18, pad=20)
ax_clo.set_xlabel("")
ax_clo.set_ylabel("")
ax_clo.tick_params(left=False, bottom=True, labelbottom=True)
plt.setp(ax_clo.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor', fontsize=18)
ax_clo.xaxis.set_tick_params(pad=5)
ax_clo.tick_params(axis='x', which='both', bottom=True, length=10, width=1.5, direction='out')

# Right panel: Subclonal heatmap + colorbar
ax_sub = plt.subplot(gs[1])
cbar_ax = plt.subplot(gs[2])
sns.heatmap(ordered_sub, cmap="Purples", ax=ax_sub,
            xticklabels=True, yticklabels=False, vmin=0, vmax=1,
            cbar=True, cbar_ax=cbar_ax)
ax_sub.set_title("Subclonal Signatures", fontsize=18, pad=20)
ax_sub.set_xlabel("")
ax_sub.set_ylabel("")
ax_sub.tick_params(left=False, bottom=True, labelbottom=True)
plt.setp(ax_sub.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor', fontsize=18)
ax_sub.xaxis.set_tick_params(pad=5)
ax_sub.tick_params(axis='x', which='both', bottom=True, length=10, width=1.5, direction='out')

# Visual divider between panels
fig.lines.append(plt.Line2D([0.5, 0.5], [0.1, 0.92], transform=fig.transFigure, color="black", linewidth=3))

# Colorbar labeling + main title and N
cbar_ax.set_title("Intensity", fontsize=18, pad=20)
cbar_ax.tick_params(labelsize=16)
fig.suptitle("Clonal and Subclonal SBS Signatures extracted with SigProfilerAssignment", fontsize=21, y=1.08)
fig.text(0.5, 1.01, "N = 277", ha='center', va='top', fontsize=21)

plt.tight_layout()
plt.savefig("ClonalvsSubclonal_SigAssign3.png", dpi=300, bbox_inches='tight')
# -------------------------------------------------------------------------------


# ==================== Figure 2: Subclonal − Clonal delta (Wilcoxon) ====================
# Reload (kept as in original script)
df_clo = pd.read_csv(
    "C:/Users/hesam/Desktop/irp5/277/Sig/SigAssign_3.4_2/Clonal/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",
    sep="\t", index_col=0
)
df_sub = pd.read_csv(
    "C:/Users/hesam/Desktop/irp5/277/Sig/SigAssign_3.4_2/Subclonal/Assignment_Solution/Activities/Assignment_Solution_Activities.txt",
    sep="\t", index_col=0
)

# Signatures to compare (fixed set)
signatures_of_interest = ['SBS1','SBS2','SBS3','SBS5','SBS6','SBS15','SBS16','SBS19','SBS30','SBS32','SBS39','SBS42','SBS44','SBS54','SBS84']

# Align samples across panels
df_clo = df_clo.sort_index()
df_sub = df_sub.sort_index()
common_samples = df_clo.index.intersection(df_sub.index)
df_clo = df_clo.loc[common_samples]
df_sub = df_sub.loc[common_samples]

# Row-wise max normalization (all columns), then select signatures of interest
df_clo_norm_all = df_clo.div(df_clo.max(axis=1), axis=0).fillna(0)
df_sub_norm_all = df_sub.div(df_sub.max(axis=1), axis=0).fillna(0)
df_clo_norm = df_clo_norm_all[signatures_of_interest]
df_sub_norm = df_sub_norm_all[signatures_of_interest]

# Mean delta per signature (Subclonal − Clonal)
signature_deltas = (df_sub_norm.mean() - df_clo_norm.mean()).sort_values(ascending=False)

# Paired Wilcoxon per signature
pvals = {}
for sig in signature_deltas.index:
    try:
        stat, p = wilcoxon(df_sub_norm[sig], df_clo_norm[sig], alternative='two-sided')
        pvals[sig] = p
    except ValueError:
        pvals[sig] = float('nan')

# Compact scientific notation for significant P values (LaTeX style)
pval_lines = []
for sig in signature_deltas.index:
    p = pvals[sig]
    if pd.notna(p) and p < 0.05:
        sci_str = f"{p:.1e}".replace("e-0", "e-").replace("e-00", "e-")
        base, exp = sci_str.split('e')
        line = rf"$P = {base}\times10^{{{int(exp)}}}$ ({sig})"
        pval_lines.append(line)
pval_text = "\n".join(pval_lines)

# Bar plot: purple = higher in Subclonal; blue = higher in Clonal
fig, ax = plt.subplots(figsize=(13, 6.2))
colors = ['purple' if v > 0 else 'blue' for v in signature_deltas]
signature_deltas.plot(kind='bar', color=colors, ax=ax)
ax.axhline(0, color='black', linewidth=0.8)
ax.set_ylabel("Subclonal - Clonal Mean Intensity", fontsize=18)
ax.set_title("SBS Signature Activity Change (Subclonal - Clonal) extracted with SigProfilerAssignment",
             fontsize=21, y=1.1, x=0.6)

# X-tick styling
plt.setp(ax.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor', fontsize=22)
ax.xaxis.set_tick_params(pad=5)
ax.tick_params(axis='x', which='both', bottom=True, length=10, width=1.5, direction='out')

# Legend + N box + significant P box
green_patch = mpatches.Patch(color='purple', label='Higher in Subclonal')
red_patch = mpatches.Patch(color='blue', label='Higher in Clonal')
ax.legend(handles=[green_patch, red_patch], loc='lower right', bbox_to_anchor=(0.4, 0.027),
          fontsize=15, frameon=False, ncol=1)
ax.text(0.98, 0.95, "N = 277", transform=ax.transAxes, fontsize=17,
        va='top', ha='right', bbox=dict(boxstyle="square,pad=0.3", edgecolor='black', facecolor='white'))
ax.text(1.02, 0.925, "Significant Differences\n" + pval_text, transform=ax.transAxes, fontsize=15,
        va='top', ha='left', bbox=dict(boxstyle="round,pad=0.5", facecolor='lightyellow', edgecolor='gray'))
ax.tick_params(axis='y', labelsize=16)

fig.subplots_adjust(left=0.145, bottom=0.25, right=0.74, top=0.85, wspace=0.2, hspace=0.2)
fig.savefig("selected_signature_delta_wilcoxon_clean_latex_FIXED_NORMALIZATION3.png", dpi=300)
# -------------------------------------------------------------------------------


# ========================================
# Output directory
outdir = "SigAssignment_SBS_Figures"
os.makedirs(outdir, exist_ok=True)

# Apply SAME style of column filtering as specified (thresholds unchanged)
keep_cols_clo = df_clo.sum(axis=0) > 100
keep_cols_sub = df_sub.sum(axis=0) > 100
df_clo_f = df_clo.loc[:, keep_cols_clo]
df_sub_f = df_sub.loc[:, keep_cols_sub]

# Use signatures present in BOTH panels
SIGS = sorted(set(df_clo_f.columns) & set(df_sub_f.columns))
print("Signatures used:", SIGS)

# Proportional (row-sum) normalization AFTER filtering
df_clo_prop = df_clo_f.div(df_clo_f.sum(axis=1), axis=0).fillna(0)
df_sub_prop = df_sub_f.div(df_sub_f.sum(axis=1), axis=0).fillna(0)

# --------------- Figure 3: Prevalence (>5%) barplot ---------------
thr = 0.05
prev_clo = (df_clo_prop.reindex(columns=SIGS, fill_value=0) > thr).mean().sort_index() * 100
prev_sub = (df_sub_prop.reindex(columns=SIGS, fill_value=0) > thr).mean().sort_index() * 100

fig, ax = plt.subplots(figsize=(13, 6.2))
x = range(len(SIGS))
ax.bar(x, prev_clo, width=0.4, label='Clonal', color='navy')
ax.bar([i+0.4 for i in x], prev_sub, width=0.4, label='Subclonal', color='purple')

ax.set_xticks([i+0.2 for i in x])
ax.set_xticklabels(SIGS)
plt.setp(ax.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor', fontsize=18)
ax.set_ylabel("Prevalence (% of samples >5%)", fontsize=18, labelpad=20)
ax.set_title("SBS Signature Prevalence (Clonal vs Subclonal) Extracted with SigProfilerAssignment",
             fontsize=18, y=1.08)
ax.legend(frameon=False, fontsize=15, loc='upper right')
ax.tick_params(axis='y', labelsize=16)
ax.text(0.015, 0.96, "N = 277", transform=ax.transAxes, fontsize=14,
        va='top', bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", facecolor="white"))

fig.subplots_adjust(bottom=0.25)
fig.savefig(os.path.join(outdir, "SBS_prevalence_barplot.png"), dpi=300)
plt.show()

# --------------- Figure 4: Mean contribution barplot ---------------
mean_clo = df_clo_prop.reindex(columns=SIGS, fill_value=0).mean().sort_index()
mean_sub = df_sub_prop.reindex(columns=SIGS, fill_value=0).mean().sort_index()

fig, ax = plt.subplots(figsize=(13, 6.2))
x = range(len(SIGS))
ax.bar(x, mean_clo, width=0.4, label='Clonal', color='steelblue')
ax.bar([i+0.4 for i in x], mean_sub, width=0.4, label='Subclonal', color='mediumpurple')
ax.axhline(0, color='black', linewidth=0.8)

ax.set_xticks([i+0.2 for i in x])
ax.set_xticklabels(SIGS)
plt.setp(ax.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor', fontsize=18)
ax.set_ylabel("Mean Contribution (proportion)", fontsize=18, labelpad=20)
ax.set_title("Mean SBS Signature Contribution (Clonal vs Subclonal) Extracted with SigProfilerAssignment",
             fontsize=18, y=1.08)
ax.legend(frameon=False, fontsize=15, loc='upper right')
ax.tick_params(axis='y', labelsize=16)
ax.text(0.015, 0.96, "N = 277", transform=ax.transAxes, fontsize=14,
        va='top', bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", facecolor="white"))

fig.subplots_adjust(bottom=0.25)
fig.savefig(os.path.join(outdir, "SBS_mean_contribution_barplot.png"), dpi=300)
plt.show()

# --------------- Figure 5: Co-occurrence (Spearman) heatmaps ---------------
corr_clo = df_clo_prop.reindex(columns=SIGS, fill_value=0).corr(method='spearman')
corr_sub = df_sub_prop.reindex(columns=SIGS, fill_value=0).corr(method='spearman')

fig = plt.figure(figsize=(20, 9))
gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 0.04], wspace=0.3)
ax1 = fig.add_subplot(gs[0])
ax2 = fig.add_subplot(gs[1])
cax = fig.add_subplot(gs[2])

sns.heatmap(corr_clo, vmin=-1, vmax=1, cmap="coolwarm", square=True,
            cbar=False, ax=ax1, xticklabels=True, yticklabels=True)
plt.setp(ax1.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor', fontsize=19)
plt.setp(ax1.get_yticklabels(), fontsize=19)
ax1.set_title("Clonal", fontsize=21, pad=12)

sns.heatmap(corr_sub, vmin=-1, vmax=1, cmap="coolwarm", square=True,
            cbar=False, ax=ax2, xticklabels=True, yticklabels=True)
plt.setp(ax2.get_xticklabels(), rotation=45, ha='right', rotation_mode='anchor', fontsize=19)
plt.setp(ax2.get_yticklabels(), fontsize=19)
ax2.set_title("Subclonal", fontsize=21, pad=12)

norm = mpl.colors.Normalize(vmin=-1, vmax=1)
sm = mpl.cm.ScalarMappable(cmap="coolwarm", norm=norm)
sm.set_array([])
cb = fig.colorbar(sm, cax=cax)
cb.set_label("Intensity", fontsize=21, labelpad=14)
cb.ax.tick_params(labelsize=18)

fig.suptitle("SBS Signature Co-occurrence (Spearman Correlation) Extracted with SigProfilerAssignment (N=277)",
             fontsize=26, y=0.97)

plt.tight_layout(rect=[0, 0, 1, 0.96])
plt.savefig(os.path.join(outdir, "SBS_correlation_heatmaps_large.png"), dpi=300)
plt.show()

# Significant pairwise Spearman (FDR < 0.05) — save full and significant tables
def bh_adjust(p):
    p = np.asarray(p)
    n = p.size
    order = np.argsort(p)
    ranks = np.arange(1, n + 1)
    q_sorted = (p[order] * n / ranks)
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q = np.empty_like(q_sorted)
    q[order] = q_sorted
    return q

def sig_spearman(df, cols, alpha=0.05):
    rows = []
    for a, b in itertools.combinations(cols, 2):
        rho, p = spearmanr(df[a], df[b], nan_policy='omit')
        rows.append((a, b, rho, p))
    res = pd.DataFrame(rows, columns=["sig1", "sig2", "rho", "p"])
    res["q"] = bh_adjust(res["p"].values)
    res = res.sort_values("q")
    return res[res["q"] < alpha], res

clo_mat = df_clo_prop.reindex(columns=SIGS, fill_value=0)
sub_mat = df_sub_prop.reindex(columns=SIGS, fill_value=0)
sig_clo, all_clo = sig_spearman(clo_mat, SIGS, alpha=0.05)
sig_sub, all_sub = sig_spearman(sub_mat, SIGS, alpha=0.05)

print("\n== Significant Spearman correlations (Clonal, FDR<0.05) ==")
print(sig_clo.head(15).to_string(index=False, formatters={"rho": "{:.3f}".format, "q": "{:.2e}".format}))
print("\n== Significant Spearman correlations (Subclonal, FDR<0.05) ==")
print(sig_sub.head(15).to_string(index=False, formatters={"rho": "{:.3f}".format, "q": "{:.2e}".format}))

all_clo.to_csv(os.path.join(outdir, "Spearman_all_pairs_clonal.csv"), index=False)
all_sub.to_csv(os.path.join(outdir, "Spearman_all_pairs_subclonal.csv"), index=False)
sig_clo.to_csv(os.path.join(outdir, "Spearman_significant_clonal_FDR05.csv"), index=False)
sig_sub.to_csv(os.path.join(outdir, "Spearman_significant_subclonal_FDR05.csv"), index=False)
# -------------------------------------------------------------------------------
