from joblib import Parallel, delayed
import os
import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy.special
import scipy.spatial
import scipy.stats
import scipy.linalg
import scipy.optimize
import scipy.integrate
from clonesig.run_clonesig import get_MU, run_clonesig
import scipy as sp
import seaborn as sns
from seaborn.distributions import _freedman_diaconis_bins
import matplotlib.ticker as ticker

#Global setup
MU = get_MU(cosmic_version=3, cancer_type=19, exome=True, artefact=False)
contexts = [str(i) for i in range(MU.shape[1])]
master_output_dir = f"clonesig_analysis_{time.strftime('%Y%m%d_%H%M%S')}"
os.makedirs(master_output_dir, exist_ok=True)

#Load input table
try:
    df_full = pd.read_csv("table9.tsv", sep="\t")
except FileNotFoundError:
    print("Error: table9.tsv not found.")
    exit()

if "sample_id" not in df_full.columns:
    print("Error: 'sample_id' column not found.")
    exit()

unique_sample_ids = df_full["sample_id"].unique()

#Per-sample CloneSig wrapper (for parallel run)
def run_clonesig_for_sample(sample_id):
    print(f"--- Processing Sample: {sample_id} ---")
    sample_output_dir = os.path.join(master_output_dir, sample_id)
    os.makedirs(sample_output_dir, exist_ok=True)

    # Prepare CloneSig inputs
    df_sample = df_full[df_full["sample_id"] == sample_id].copy()
    T = df_sample["context_index"].values
    B = df_sample["t_alt_count"].values
    D = df_sample["t_ref_count"].values + B
    C_normal = np.full_like(D, 2)
    C_tumor_tot = df_sample["cnTotal"].values
    C_tumor_minor = df_sample["cnMinor"].values
    purity = float(df_sample["purity"].iloc[0])

    # Run CloneSig
    try:
        est, lr, pval, new_inputMU, cst_est, future_sigs = run_clonesig(
            T=T, B=B, D=D,
            C_normal=C_normal, C_tumor_tot=C_tumor_tot, C_tumor_minor=C_tumor_minor,
            purity=purity, inputMU=MU, max_nb_clones=6,
            nb_fits=20, prefit_signatures=True, prefit_thresh=0.05
        )
    except Exception as e:
        print(f"Error running CloneSig for {sample_id}: {e}")
        return

    # Extract estimates and assignments
    qun, vmnu, rnus = est.qun, est.vmnu, est.rnus
    clone_assignments = np.argmax(qun, axis=1)

    # Map to COSMIC signature indices
    if isinstance(new_inputMU, pd.DataFrame):
        signature_indices = list(map(int, new_inputMU.columns))
    else:
        signature_indices = list(range(new_inputMU.shape[0]))

    # Assemble per-mutation table
    est_table = pd.DataFrame({
        "sample_id": sample_id,
        "trinucleotide": [contexts[i] for i in est.T],
        "var_counts": est.B,
        "minor_cn": est.C_tumor_minor,
        "major_cn": est.C_tumor_major,
        "total_cn": est.C_tumor_tot,
        "depth": est.D,
        "clonesig_clone_reassignment": clone_assignments,
        "signature": [signature_indices[np.argmax(rnus[i, clone_assignments[i], :])] for i in range(est.N)],
        "mult": [np.argmax(vmnu[i, clone_assignments[i], :]) + 1 for i in range(est.N)]
    })

    # Compute VAF-derived metrics
    vaf = est_table["var_counts"].values / est_table["depth"].values
    total_cn = est_table["total_cn"].values
    mult = est_table["mult"].values
    est_table["vaf"] = vaf
    est_table["vaf_cn"] = vaf * total_cn / mult
    est_table["vaf_purity"] = vaf / est.p * (((1 - est.p) * 2 + est.p * total_cn) / mult)

    est_table.to_csv(
        os.path.join(sample_output_dir, f"{sample_id}_clonesig_est_table.csv"),
        index=False
    )

    # Plot: purity-adjusted VAF distribution by clone
    sns.set_context("talk")
    nb_bins = min(_freedman_diaconis_bins(est_table.vaf_purity) * 2, 50)
    plt.figure(figsize=(10, 6))
    sns.histplot(est_table.vaf_purity, bins=nb_bins, label="All Mutations", kde=False, color="lightgray")
    for i in range(est.J):
        sns.histplot(
            est_table[est_table.clonesig_clone_reassignment == i].vaf_purity,
            bins=nb_bins, kde=False, label=f"Clone {i}", alpha=0.6
        )
    plt.title(f"VAF Purity Distribution – {sample_id}")
    plt.xlabel("VAF (Purity-adjusted)")
    plt.ylabel("Number of Mutations")
    plt.legend(title="Clones")
    plt.tight_layout()
    plt.savefig(os.path.join(sample_output_dir, f"{sample_id}_vaf_purity_clonesig_reconstruction.png"))
    plt.close()

    # Plot: signature activity per clone
    plt.figure(figsize=(12, 3))
    ax = sns.heatmap(
        est.pi,
        vmin=0, vmax=1, cmap="GnBu",
        yticklabels=[f"Clone {i}" for i in range(est.J)],
        xticklabels=signature_indices
    )
    ax.set_ylabel("Clone index", fontsize=12)
    ax.set_xlabel("COSMIC Signature index (0–64)", fontsize=12)
    ax.set_title(f"{sample_id} – signature activity per clone (CloneSig)", fontsize=14)
    for idx, label in enumerate(ax.get_xticklabels()):
        if idx % 2 != 0:
            label.set_visible(False)
        else:
            label.set_fontsize(8)
            label.set_rotation(90)
    plt.tight_layout()
    plt.savefig(os.path.join(sample_output_dir, f"{sample_id}_clonesig_signature_activity_per_clone.png"), dpi=300)
    plt.close()

#Parallel execution
Parallel(n_jobs=4)(delayed(run_clonesig_for_sample)(sid) for sid in unique_sample_ids)
print(f"\nAll {len(unique_sample_ids)} samples processed in parallel. Outputs in: {master_output_dir}")
