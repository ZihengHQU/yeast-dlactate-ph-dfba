"""Assemble the source-data package behind every figure."""

import json
import os
import shutil

import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SRC = os.path.join(ROOT, "source_data")
os.makedirs(SRC, exist_ok=True)

PKA = 3.86


def fig2_curves():
    p = json.load(open(os.path.join(ROOT, "data", "toxicity_params.json")))
    m50, k_h, k_a = p["m50"], p["k_h"], p["k_a"]
    ha5 = 900 * (1 / (1 + 10 ** (5.0 - PKA)))
    k_pure = m50 / ha5
    ph = np.round(np.arange(2.6, 5.61, 0.05), 2)
    f = 1 / (1 + 10 ** (ph - PKA))
    rows = pd.DataFrame({
        "pH": ph,
        "total_acid_for_50pct_inhibition_single_mM": m50 / (k_pure * f),
        "total_acid_for_50pct_inhibition_two_mM": m50 / (k_h * f + k_a * (1 - f)),
    })
    rows.to_csv(os.path.join(SRC, "figure2_toxicity_curves.csv"), index=False)

    obs = pd.DataFrame({
        "pH": [5.0, 3.5, 3.0],
        "total_acid_mM": [900, 750, 500],
        "role": ["calibration", "held-out test", "calibration"],
    })
    obs["undissociated_mM"] = obs["total_acid_mM"] * (1 / (1 + 10 ** (obs["pH"] - PKA)))
    obs["anion_mM"] = obs["total_acid_mM"] - obs["undissociated_mM"]
    obs.to_csv(os.path.join(SRC, "figure2_observations.csv"), index=False)


def copy(src_name, dst_name, note):
    src = os.path.join(ROOT, "results", src_name)
    if os.path.exists(src):
        shutil.copy(src, os.path.join(SRC, dst_name))
        return dst_name, note
    return None, note


def main():
    fig2_curves()

    manifest = [
        ("figure1_mechanism", "-", "Conceptual diagram; no numerical data."),
        ("figure2_toxicity", "figure2_observations.csv",
         "Abbott et al. 2008 chemostat points (calibration and held-out test)."),
        ("figure2_toxicity", "figure2_toxicity_curves.csv",
         "Model curves plotted in Fig. 2A (recomputed from data/toxicity_params.json)."),
        ("figure3_batch", "figure3_traj_Baek_batch.csv", "Uncontrolled-pH batch trajectory."),
        ("figure3_batch", "figure3_traj_Ishida_batch.csv", "pH 6.5 batch trajectory."),
        ("figure4_validation", "figure4_counterfactual.csv",
         "Predictions of both toxicity models against eight strains."),
        ("figure5_design", "figure5_ph_sweep.csv", "Titre versus controlled pH."),
        ("figure5_design", "figure5_genotypes.csv", "Titre by genotype at two pH values."),
        ("figure6_sensitivity", "figure6_sensitivity.csv",
         "Titre versus LDH capacity and oxygen supply."),
    ]

    copy("traj_Baek_batch.csv", "figure3_traj_Baek_batch.csv", "")
    copy("traj_Ishida_batch.csv", "figure3_traj_Ishida_batch.csv", "")
    copy("counterfactual.csv", "figure4_counterfactual.csv", "")
    copy("ph_sweep.csv", "figure5_ph_sweep.csv", "")
    copy("genotypes.csv", "figure5_genotypes.csv", "")
    copy("sensitivity.csv", "figure6_sensitivity.csv", "")

    pd.DataFrame(manifest, columns=["figure", "file", "description"]).to_csv(
        os.path.join(SRC, "MANIFEST.csv"), index=False)
    print("source data files:")
    for f in sorted(os.listdir(SRC)):
        print("  ", f, os.path.getsize(os.path.join(SRC, f)))


if __name__ == "__main__":
    main()
