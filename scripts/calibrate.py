"""Calibrate the two-term lactic-acid toxicity module.

Data: Abbott, Suir, van Maris & Pronk (2008) AEM 74:5759-5768.
Anaerobic glucose-limited chemostat, D = 0.10 h^-1, CEN.PK 113-7D.
Total lactic acid that halves the biomass yield: 900 mM (pH 5.0), 500 mM (pH 3.0);
the pH 3.5 point (750 mM) is held out as a test of the fitted model.
"""

import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from src import core

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
os.makedirs(os.path.join(ROOT, "data"), exist_ok=True)

FIT = [(5.0, 900.0), (3.0, 500.0)]
TEST = (3.5, 750.0)


def m50_for_half_yield(model, D=0.1):
    y0 = core.chemostat_yield(model, D=D, maintenance=0.0)
    lo, hi = 0.0, 200.0
    for _ in range(45):
        mid = 0.5 * (lo + hi)
        y = core.chemostat_yield(model, D=D, maintenance=mid)
        if np.isnan(y) or y < 0.5 * y0:
            hi = mid
        else:
            lo = mid
    return 0.5 * (lo + hi), y0


def main():
    m = core.load_model()
    m50, y0 = m50_for_half_yield(m)
    print(f"Y0(D=0.1) = {y0:.4f} g/g ; m50 = {m50:.2f} mmol ATP gDW-1 h-1")

    A, rhs = [], []
    for pH, C in FIT:
        ha, an = core.speciation(C, pH)
        A.append([ha, an])
        rhs.append(m50)
    k_h, k_a = np.linalg.solve(np.array(A), np.array(rhs))
    print(f"k_h = {k_h:.5f}  k_a = {k_a:.6f} (mmol ATP per mmol)")

    for pH, C in FIT + [TEST]:
        ha, an = core.speciation(C, pH)
        pred = k_h * ha + k_a * an
        print(f"  pH={pH:4.1f} C={C:4.0f}  HA={ha:6.1f} A-={an:6.1f} -> m={pred:6.2f}")

    # concentration that halves yield at pH 3.5 according to the model
    ha_u, an_u = core.speciation(1.0, TEST[0])
    c_pred = m50 / (k_h * ha_u + k_a * an_u)
    print(f"pH 3.5: model predicts 50%% effect at {c_pred:.0f} mM, observed {TEST[1]:.0f} mM "
          f"({(c_pred / TEST[1] - 1) * 100:+.0f}%)")
    # what a pure undissociated-acid model would predict at pH 3.5
    c_ha = m50 / (k_h * ha_u)
    print(f"undissociated-only model would predict {c_ha:.0f} mM at pH 3.5 "
          f"(observed {TEST[1]:.0f} mM, i.e. {TEST[1] / c_ha:.1f}x higher)")

    out = dict(k_h=k_h, k_a=k_a, m50=m50, Y0=y0,
               fit=[dict(pH=p, C=c) for p, c in FIT],
               test=dict(pH=TEST[0], C=TEST[1], C_model=c_pred))
    with open(os.path.join(ROOT, "data", "toxicity_params.json"), "w") as fh:
        json.dump(out, fh, indent=2)
    print("written data/toxicity_params.json")


if __name__ == "__main__":
    main()
