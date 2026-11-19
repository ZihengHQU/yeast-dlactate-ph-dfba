"""Expanded validation (8 strains) and the single-mechanism counterfactual.

Two toxicity models are applied to the same strain set with no parameter fitted
to any titre:
  two-mechanism   k_H, k_A calibrated on the pH 5.0 and pH 3.0 chemostats
  single-mechanism k_H only, calibrated on the pH 5.0 chemostat
"""

import json
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src import core

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
os.makedirs(RES, exist_ok=True)

Q_O2, Q_LAC, BUF = 5.0, 8.0, 175.0

STRAINS = [
    dict(name="Ishida 2006", del_=["PDC1"], copies=2.0, mode="batch", pH=6.5,
         s0=100.0, t_end=140, titer=61.5),
    dict(name="Rahmasari 2025 (PDC1)", del_=["PDC1"], copies=1.0, mode="batch", pH=None,
         s0=100.0, t_end=140, titer=51.0),
    dict(name="Rahmasari 2025 (CYB2)", del_=["PDC1"], copies=1.0, mode="batch", pH=None,
         s0=100.0, t_end=140, titer=31.9),
    dict(name="Baek 2016 (no pH control)", del_=["DLD1", "JEN1", "PDC1", "ADH1", "GPD1", "GPD2"],
         copies=1.0, mode="batch", pH=None, s0=100.0, t_end=140, titer=48.9),
    dict(name="Baek 2016 (pH 6)", del_=["DLD1", "JEN1", "PDC1", "ADH1", "GPD1", "GPD2"],
         copies=2.0, mode="fedbatch", pH=6.0, s0=30.0, t_end=140, titer=112.0),
    dict(name="Baek 2017 (pH 3.5)", del_=["DLD1", "PDC1", "ADH1", "GPD1", "GPD2"],
         copies=2.0, mode="fedbatch", pH=3.5, s0=30.0, t_end=140, titer=82.6),
    dict(name="Zhong 2019 (no pH control)", del_=["PDC1", "PDC6", "JEN1", "CYB2", "DLD1", "ADH1"],
         copies=3.0, mode="fedbatch", pH=None, s0=30.0, t_end=160, titer=80.0),
    dict(name="Sornlek 2022 hybrid2", del_=["PDC1", "CYB2", "GPD1", "ADH1", "ALD6"],
         copies=1.0, mode="batch", pH=None, s0=100.0, t_end=140, titer=23.4),
]


def make(case):
    m = core.load_model()
    core.apply_genotype(m, case["del_"], case["copies"])
    m.reactions.r_1992_REV.bounds = (0.0, Q_O2)
    m.reactions.r_1546.upper_bound = Q_LAC
    return m


def run(case, tox):
    m = make(case)
    if case["mode"] == "batch":
        return core.simulate_batch(m, pH=case["pH"] or 6.0, X0=0.5,
                                   S0=case["s0"] / core.MW_GLUCOSE, t_end=case["t_end"],
                                   dt=0.25, tox=tox,
                                   drift_buffer=None if case["pH"] else BUF)
    return core.simulate_fedbatch(m, pH=case["pH"] or 6.0, X0=0.5,
                                  S0=case["s0"] / core.MW_GLUCOSE, V0=0.5,
                                  t_end=case["t_end"], dt=0.25, tox=tox,
                                  S_feed=550.0, S_target=20.0,
                                  drift_buffer=None if case["pH"] else BUF)


def main():
    p = json.load(open(os.path.join(ROOT, "data", "toxicity_params.json")))
    m50 = p["m50"]
    ha5, _ = core.speciation(900.0, 5.0)
    tox_two = dict(k_h=p["k_h"], pK_phi=None, k_a=p["k_a"])
    tox_one = dict(k_h=m50 / ha5, pK_phi=None, k_a=0.0)

    rows = []
    for case in STRAINS:
        for tag, tox in (("two", tox_two), ("single", tox_one)):
            core.TOX.clear(); core.TOX.update(tox)
            df = run(case, tox)
            rows.append(dict(strain=case["name"], model=tag,
                             titer=df["L"].iloc[-1] * core.MW_LACTATE,
                             titer_lit=case["titer"], time_h=df["t"].iloc[-1]))
            print(f"{case['name']:28s} {tag:6s} {rows[-1]['titer']:6.1f} "
                  f"(lit {case['titer']:6.1f})", flush=True)
    out = pd.DataFrame(rows)
    out["err_pct"] = (out["titer"] - out["titer_lit"]) / out["titer_lit"] * 100
    out.to_csv(os.path.join(RES, "counterfactual.csv"), index=False)

    for tag in ("two", "single"):
        sub = out[out["model"] == tag]
        print(f"{tag}: MAE = {np.mean(np.abs(sub.err_pct)):.1f}%  bias = {np.mean(sub.err_pct):+.1f}%")


if __name__ == "__main__":
    main()
