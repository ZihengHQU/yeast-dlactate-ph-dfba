"""Full simulation suite: validation, pH sweeps, genotype comparison, sensitivity.

Everything is driven by data/toxicity_params.json (calibrated on chemostat physiology).
O2 supply and LDH capacity are physical constraints, not fitted to titres:
  q_O2_max   = 5 mmol gDW-1 h-1  (kLa-limited upper range for lab bioreactors)
  q_lac_max  = 8 mmol gDW-1 h-1  (~0.5% of proteome at kcat 144 h-1)
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

Q_O2 = 5.0
Q_LAC = 8.0
BUF = 175.0


def tox():
    with open(os.path.join(ROOT, "data", "toxicity_params.json")) as fh:
        p = json.load(fh)
    core.TOX.update(k_h=p["k_h"], pK_phi=None, k_a=p["k_a"])
    return core.TOX


def make(deletions, copies=1.0, q_o2=Q_O2, q_lac=Q_LAC):
    m = core.load_model()
    core.apply_genotype(m, deletions, copies)
    m.reactions.r_1992_REV.bounds = (0.0, q_o2)
    if q_lac is not None:
        m.reactions.r_1546.upper_bound = q_lac
    return m


STRAINS = [
    dict(name="Ishida 2006", del_=["PDC1"], copies=2.0, mode="batch", pH=6.5,
         s0=100.0, t_end=120, titer=61.5, yld=0.612),
    dict(name="Rahmasari 2025", del_=["PDC1"], copies=1.0, mode="batch", pH=None,
         s0=100.0, t_end=120, titer=51.0, yld=np.nan),
    dict(name="Baek 2016 (non-neutralised)", del_=["DLD1", "JEN1", "PDC1", "ADH1", "GPD1", "GPD2"],
         copies=1.0, mode="batch", pH=None, s0=100.0, t_end=120, titer=48.9, yld=np.nan),
    dict(name="Baek 2016 (neutralised)", del_=["DLD1", "JEN1", "PDC1", "ADH1", "GPD1", "GPD2"],
         copies=2.0, mode="fedbatch", pH=6.0, s0=30.0, t_end=140, titer=112.0, yld=0.80),
    dict(name="Baek 2017 (pH 3.5)", del_=["DLD1", "PDC1", "ADH1", "GPD1", "GPD2"],
         copies=2.0, mode="fedbatch", pH=3.5, s0=30.0, t_end=140, titer=82.6, yld=0.83),
]


def run(s, q_lac=Q_LAC, q_o2=Q_O2, store=None):
    m = make(s["del_"], s["copies"], q_o2=q_o2, q_lac=q_lac)
    if s["mode"] == "batch":
        df = core.simulate_batch(m, pH=s["pH"] or 6.0, X0=0.5,
                                 S0=s["s0"] / core.MW_GLUCOSE, t_end=s["t_end"], dt=0.25,
                                 tox=core.TOX, drift_buffer=None if s["pH"] else BUF)
    else:
        df = core.simulate_fedbatch(m, pH=s["pH"], X0=0.5, S0=s["s0"] / core.MW_GLUCOSE,
                                    V0=0.5, t_end=s["t_end"], dt=0.25, tox=core.TOX,
                                    S_feed=550.0, S_target=20.0)
    if store is not None:
        df.to_csv(store, index=False)
    return df


def main():
    tox()
    rows = []
    for s in STRAINS:
        df = run(s, store=os.path.join(RES, f"traj_{s['name'].split()[0]}_{s['mode']}.csv"))
        titer = df["L"].iloc[-1] * core.MW_LACTATE
        glc0 = s["s0"] / core.MW_GLUCOSE
        used = (glc0 - (df["S"].iloc[-1] if s["mode"] == "batch" else 0.0)) * core.MW_GLUCOSE
        rows.append(dict(strain=s["name"], mode=s["mode"], pH=s["pH"],
                         titer_model=titer, titer_lit=s["titer"],
                         err_pct=(titer - s["titer"]) / s["titer"] * 100,
                         yield_model=titer / used if s["mode"] == "batch" and used > 1 else np.nan,
                         yield_lit=s["yld"],
                         time_h=df["t"].iloc[-1],
                         prod_model=titer / df["t"].iloc[-1]))
        print(f"{s['name']:32s} {titer:7.1f} g/L  ({rows[-1]['err_pct']:+.0f}%)  lit {s['titer']}")
    val = pd.DataFrame(rows)
    val.to_csv(os.path.join(RES, "validation_strains.csv"), index=False)

    # ---- pH sweep for the production strain
    sweep = []
    base = dict(del_=["DLD1", "JEN1", "PDC1", "ADH1", "GPD1", "GPD2"], copies=2.0)
    for pH in [2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0]:
        s = dict(**base, name="sweep", mode="batch", pH=pH, s0=100.0, t_end=140)
        df = run(s, store=os.path.join(RES, f"traj_pH_{pH}.csv"))
        sweep.append(dict(pH=pH, titer=df["L"].iloc[-1] * core.MW_LACTATE,
                          time=df["t"].iloc[-1], mu_end=df["mu"].iloc[-1],
                          m_acid=df["m_acid"].iloc[-1]))
        print(f"pH {pH}: {sweep[-1]['titer']:.1f} g/L")
    pd.DataFrame(sweep).to_csv(os.path.join(RES, "ph_sweep.csv"), index=False)

    # ---- genotype comparison
    genotypes = {
        "pdc1Δ": ["PDC1"],
        "pdc1Δ + gpd1/2Δ": ["PDC1", "GPD1", "GPD2"],
        "pdc1Δ + gpdΔ + dld1Δ + jen1Δ": ["PDC1", "GPD1", "GPD2", "DLD1", "JEN1"],
        "full (Baek 2016)": ["DLD1", "JEN1", "PDC1", "ADH1", "GPD1", "GPD2"],
    }
    g_rows = []
    for name, dl in genotypes.items():
        for pH in [6.0, 3.5]:
            s = dict(name=name, del_=dl, copies=2.0, mode="batch", pH=pH, s0=100.0, t_end=140)
            df = run(s)
            g_rows.append(dict(genotype=name, pH=pH,
                               titer=df["L"].iloc[-1] * core.MW_LACTATE,
                               biomass=df["X"].iloc[-1], time=df["t"].iloc[-1]))
    pd.DataFrame(g_rows).to_csv(os.path.join(RES, "genotypes.csv"), index=False)
    print(pd.DataFrame(g_rows).round(1).to_string(index=False))

    # ---- sensitivity: q_lac and q_O2
    sens = []
    for ql in [3, 5, 8, 12, 20]:
        for pH in [6.0, 3.5]:
            s = dict(**base, name="sens", mode="batch", pH=pH, s0=100.0, t_end=140)
            df = run(s, q_lac=ql)
            sens.append(dict(param="q_lac_max", value=ql, pH=pH,
                             titer=df["L"].iloc[-1] * core.MW_LACTATE))
    for qo in [2, 3, 5, 8, 12]:
        for pH in [6.0, 3.5]:
            s = dict(**base, name="sens", mode="batch", pH=pH, s0=100.0, t_end=140)
            df = run(s, q_o2=qo)
            sens.append(dict(param="q_O2_max", value=qo, pH=pH,
                             titer=df["L"].iloc[-1] * core.MW_LACTATE))
    pd.DataFrame(sens).to_csv(os.path.join(RES, "sensitivity.csv"), index=False)
    print("sensitivity done")


if __name__ == "__main__":
    main()
