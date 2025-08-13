import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import gridspec
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import matplotlib.patches as mpatches
from scipy.stats import wilcoxon, spearmanr
import itertools

# --------------------------- Paths & I/O ---------------------------
subclonal_csv = r"C:/Users/hesam/Desktop/irp5/277/Sig_CNV/Subclonal/Fit_result_Subclonal_CNV.csv"
clonal_csv    = r"C:/Users/hesam/Desktop/irp5/277/Sig_CNV/Clonal/Fit_result_Clonal_CNV.csv"
outdir        = "SigCNVFigures"
os.makedirs(outdir, exist_ok=True)
sns.set(context="talk")

# --------------------------- Load & preprocess ---------------------------
def load_clean(path, min_sum=15):
    df = pd.read_csv(path, index_col=0)
    df = df[df.sum(axis=1) >= min_sum]                  # drop weak samples
    df = df.apply(pd.to_numeric, errors="coerce").fillna(0).T  # samples x signatures
    return df

df_sub = load_clean(subclonal_csv)
df_clo = load_clean(clonal_csv)

# harmonizing samples
common_samples = df_sub.index.intersection(df_clo.index)
df_sub = df_sub.loc[common_samples]
df_clo = df_clo.loc[common_samples]
N = len(common_samples)

# dropping unwanted signatures
excluded = ['CN3','CN5','CN6','CN7','CN10','CN12','CN14','CN18']
df_sub = df_sub.drop(columns=[s for s in excluded if s in df_sub.columns], errors="ignore")
df_clo = df_clo.drop(columns=[s for s in excluded if s in df_clo.columns], errors="ignore")

# keeping common signatures and row-normalize (proportions)
SIGS = sorted(set(df_sub.columns) & set(df_clo.columns))
df_sub_prop = df_sub.reindex(columns=SIGS).div(df_sub.reindex(columns=SIGS).sum(axis=1), axis=0).fillna(0)
df_clo_prop = df_clo.reindex(columns=SIGS).div(df_clo.reindex(columns=SIGS).sum(axis=1), axis=0).fillna(0)

# --------------------------- Heatmaps: clonal vs subclonal ---------------------------
def heatmaps_clonal_subclonal(df_clo_norm, df_sub_norm, out_png):
    fig = plt.figure(figsize=(14, 6))
    gs = gridspec.GridSpec(1, 3, width_ratios=[5, 5, 0.3], wspace=0.03)

    ax_tr = plt.subplot(gs[0])
    sns.heatmap(df_clo_norm, cmap="Blues", ax=ax_tr, xticklabels=True, yticklabels=False, vmin=0, vmax=1, cbar=False)
    ax_tr.set_title("Clonal CNV Signatures", fontsize=18, pad=20)
    plt.setp(ax_tr.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor", fontsize=18)
    ax_tr.xaxis.set_tick_params(pad=5)
    ax_tr.tick_params(axis="x", which="both", bottom=True, length=10, width=1.5, direction="out")

    ax_sb = plt.subplot(gs[1])
    cbar_ax = plt.subplot(gs[2])
    sns.heatmap(df_sub_norm, cmap="Purples", ax=ax_sb, xticklabels=True, yticklabels=False, vmin=0, vmax=1,
                cbar=True, cbar_ax=cbar_ax)
    ax_sb.set_title("Subclonal CNV Signatures", fontsize=18, pad=20)
    plt.setp(ax_sb.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor", fontsize=18)
    ax_sb.xaxis.set_tick_params(pad=5)
    ax_sb.tick_params(axis="x", which="both", bottom=True, length=10, width=1.5, direction="out")

    cbar_ax.set_title("Intensity", fontsize=18, pad=20)
    cbar_ax.tick_params(labelsize=16)

    fig.lines.append(plt.Line2D([0.5, 0.5], [0.1, 0.92], transform=fig.transFigure, color="black", linewidth=3))
    fig.suptitle("Clonal and Subclonal CNV Signatures Extracted with Sigminer", fontsize=21, y=1.08)
    fig.text(0.5, 1.01, f"N = {N}", ha="center", va="top", fontsize=21)

    plt.tight_layout()
    fig.savefig(out_png, dpi=300, bbox_inches="tight")
    plt.show()

# max-normalized view
df_sub_max = df_sub.div(df_sub.max(axis=1), axis=0).fillna(0)
df_clo_max = df_clo.div(df_clo.max(axis=1), axis=0).fillna(0)
heatmaps_clonal_subclonal(df_clo_max, df_sub_max, out_png=os.path.join(outdir, "CNV_Clonal_vs_Subclonal_maxnorm.png"))

# --------------------------- Delta barplot + Wilcoxon ---------------------------
def cnv_delta_wilcoxon(df_sub_prop, df_clo_prop, sigs, out_png):
    deltas = (df_sub_prop.mean() - df_clo_prop.mean()).loc[sigs].sort_values(ascending=False)

    pvals = {}
    for sig in deltas.index:
        a = df_sub_prop[sig]
        b = df_clo_prop[sig]
        try:
            _, p = wilcoxon(a, b, alternative="two-sided")
        except ValueError:
            p = np.nan
        pvals[sig] = p
    pvals = pd.Series(pvals)

    lines = []
    for sig, p in pvals.dropna().items():
        if p < 0.05:
            s = f"{p:.1e}".replace("e-0", "e-").replace("e-00", "e-")
            base, exp = s.split("e")
            lines.append(rf"$P = {base}\times10^{{{int(exp)}}}$ ({sig})")
    p_text = "\n".join(lines)

    fig, ax = plt.subplots(figsize=(13, 6.2))
    colors = ["purple" if v > 0 else "blue" for v in deltas]
    deltas.plot(kind="bar", color=colors, ax=ax)

    ax.axhline(0, color="black", linewidth=0.8)
    ax.set_ylabel("Subclonal − Clonal Mean Intensity", fontsize=18)
    ax.set_title("CNV Signature Change (Subclonal − Clonal) extracted with Sigminer", fontsize=21, y=1.15, x=0.6)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor", fontsize=22)
    ax.tick_params(axis="x", which="both", bottom=True, length=10, width=1.5)

    purple_patch = mpatches.Patch(color="purple", label="Higher in Subclonal")
    blue_patch   = mpatches.Patch(color="blue",   label="Higher in Clonal")
    ax.legend(handles=[purple_patch, blue_patch], loc="lower right",
              bbox_to_anchor=(0.4, 0.027), fontsize=15, frameon=False, ncol=1)

    ax.text(0.98, 0.95, f"N = {N}", transform=ax.transAxes, fontsize=17,
            va="top", ha="right", bbox=dict(boxstyle="square,pad=0.3", edgecolor="black", facecolor="white"))
    ax.text(1.02, 0.925, "Significant Differences\n" + p_text, transform=ax.transAxes, fontsize=15,
            va="top", ha="left", bbox=dict(boxstyle="round,pad=0.5", facecolor="lightyellow", edgecolor="gray"))

    fig.subplots_adjust(left=0.145, bottom=0.25, right=0.74, top=0.85)
    fig.savefig(out_png, dpi=300)
    plt.show()
    return deltas, pvals

# choosing signatures to display based on first figure
signatures_of_interest = ['CN1','CN2','CN9','CN13','CN17','CN20','CN21','CN22','CN23','CN24']
signatures_of_interest = [s for s in signatures_of_interest if s in SIGS]

cnv_deltas, cnv_pvals = cnv_delta_wilcoxon(
    df_sub_prop, df_clo_prop, signatures_of_interest,
    out_png=os.path.join(outdir, "CNV_signature_delta_wilcoxon.png")
)
print("Mean CNV deltas (Subclonal − Clonal):\n", cnv_deltas)

# --------------------------- Prevalence (>5%) ---------------------------
def prevalence_barplot(df_clo_prop, df_sub_prop, sigs, thr, out_png):
    prev_clo = (df_clo_prop[sigs] > thr).mean().sort_index() * 100
    prev_sub = (df_sub_prop[sigs] > thr).mean().sort_index() * 100

    fig, ax = plt.subplots(figsize=(13, 6.2))
    x = range(len(sigs))
    ax.bar(x, prev_clo, width=0.4, label="Clonal", color="navy")
    ax.bar([i + 0.4 for i in x], prev_sub, width=0.4, label="Subclonal", color="purple")

    ax.set_xticks([i + 0.2 for i in x])
    ax.set_xticklabels(sigs)
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor", fontsize=18)
    ax.set_ylabel("Prevalence (% of samples >5%)", fontsize=18, labelpad=20)
    ax.set_title("CNV Signature Prevalence (Clonal vs Subclonal) Extracted with sigminer", fontsize=18, y=1.08)
    ax.legend(frameon=False, fontsize=15, loc="upper right")
    ax.tick_params(axis="y", labelsize=16)
    ax.text(0.72, 0.92, f"N = {N}", transform=ax.transAxes, fontsize=14,
            va="top", bbox=dict(boxstyle="round,pad=0.3", edgecolor="black", facecolor="white"))

    fig.subplots_adjust(bottom=0.25)
    fig.savefig(out_png, dpi=300)
    plt.show()

prevalence_barplot(
    df_clo_prop, df_sub_prop, SIGS, thr=0.05,
    out_png=os.path.join(outdir, "CNV_prevalence_barplot.png")
)

# --------------------------- Correlations (Spearman) ---------------------------
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

def spearman_tables(df_prop, sigs, alpha=0.05):
    rows = []
    for a, b in itertools.combinations(sigs, 2):
        rho, p = spearmanr(df_prop[a], df_prop[b], nan_policy="omit")
        rows.append((a, b, rho, p))
    res = pd.DataFrame(rows, columns=["sig1", "sig2", "rho", "p"])
    res["q"] = bh_adjust(res["p"].values)
    sig = res[res["q"] < alpha].sort_values("q")
    return sig, res

def corr_heatmaps(df_clo_prop, df_sub_prop, sigs, out_png):
    corr_clo = df_clo_prop[sigs].corr(method="spearman")
    corr_sub = df_sub_prop[sigs].corr(method="spearman")

    fig = plt.figure(figsize=(20, 9))
    gs = gridspec.GridSpec(1, 3, width_ratios=[1, 1, 0.04], wspace=0.3)

    ax1 = fig.add_subplot(gs[0])
    ax2 = fig.add_subplot(gs[1])
    cax = fig.add_subplot(gs[2])

    sns.heatmap(corr_clo, vmin=-1, vmax=1, cmap="coolwarm", square=True,
                cbar=False, ax=ax1, xticklabels=True, yticklabels=True)
    plt.setp(ax1.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor", fontsize=19)
    plt.setp(ax1.get_yticklabels(), rotation=0, ha="right", fontsize=19)
    ax1.set_title("Clonal", fontsize=21, pad=12)

    sns.heatmap(corr_sub, vmin=-1, vmax=1, cmap="coolwarm", square=True,
                cbar=False, ax=ax2, xticklabels=True, yticklabels=True)
    plt.setp(ax2.get_xticklabels(), rotation=45, ha="right", rotation_mode="anchor", fontsize=19)
    plt.setp(ax2.get_yticklabels(), rotation=0, ha="right", fontsize=19)
    ax2.set_title("Subclonal", fontsize=21, pad=12)

    norm = mpl.colors.Normalize(vmin=-1, vmax=1)
    sm = mpl.cm.ScalarMappable(cmap="coolwarm", norm=norm); sm.set_array([])
    cb = fig.colorbar(sm, cax=cax)
    cb.set_label("Intensity", fontsize=21, labelpad=14)
    cb.ax.tick_params(labelsize=18)

    fig.suptitle(f"CNV Signature Co-occurrence (Spearman) Extracted with sigminer (N={N})",
                 fontsize=26, y=0.97)
    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(out_png, dpi=300)
    plt.show()

corr_heatmaps(
    df_clo_prop, df_sub_prop, SIGS,
    out_png=os.path.join(outdir, "CNV_correlation_heatmaps.png")
)

sig_clo, all_clo = spearman_tables(df_clo_prop, SIGS, alpha=0.05)
sig_sub, all_sub = spearman_tables(df_sub_prop, SIGS, alpha=0.05)

print("\n== Significant Spearman correlations (Clonal, FDR<0.05) ==")
print(sig_clo.head(20).to_string(index=False,
      formatters={"rho":"{:.3f}".format, "q":"{:.2e}".format}))
print("\n== Significant Spearman correlations (Subclonal, FDR<0.05) ==")
print(sig_sub.head(20).to_string(index=False,
      formatters={"rho":"{:.3f}".format, "q":"{:.2e}".format}))

all_clo.to_csv(os.path.join(outdir, "CNV_Spearman_all_pairs_clonal.csv"), index=False)
all_sub.to_csv(os.path.join(outdir, "CNV_Spearman_all_pairs_subclonal.csv"), index=False)
sig_clo.to_csv(os.path.join(outdir, "CNV_Spearman_significant_clonal_FDR05.csv"), index=False)
sig_sub.to_csv(os.path.join(outdir, "CNV_Spearman_significant_subclonal_FDR05.csv"), index=False)
