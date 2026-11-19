import json, os, sys
import pandas as pd
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src import core
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
p = json.load(open(os.path.join(ROOT, "data", "toxicity_params.json")))
base = dict(k_h=p["k_h"], k_a=p["k_a"])
DL = ["DLD1", "JEN1", "PDC1", "ADH1", "GPD1", "GPD2"]
rows = []
for fh in [0.8, 1.0, 1.2]:
    for fa in [0.8, 1.0, 1.2]:
        core.TOX.update(k_h=base["k_h"] * fh, pK_phi=None, k_a=base["k_a"] * fa)
        for pH in [6.0, 3.5]:
            m = core.load_model(); core.apply_genotype(m, DL, 2.0)
            m.reactions.r_1992_REV.bounds = (0, 5.0); m.reactions.r_1546.upper_bound = 8.0
            df = core.simulate_batch(m, pH=pH, X0=0.5, S0=100 / core.MW_GLUCOSE, t_end=140,
                                     dt=0.5, tox=core.TOX)
            rows.append(dict(f_kh=fh, f_ka=fa, pH=pH, titer=df["L"].iloc[-1] * core.MW_LACTATE))
            print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(os.path.join(ROOT, "results", "param_uncertainty.csv"), index=False)
