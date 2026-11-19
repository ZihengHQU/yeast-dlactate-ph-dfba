"""Core of the pH-resolved enzyme-constrained dFBA framework.

Units
-----
fluxes        mmol gDW^-1 h^-1
concentrations mM (== mmol L^-1)
biomass        gDW L^-1
time          h
"""

from __future__ import annotations

import os
import warnings

import cobra
import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
MODEL_PATH = os.path.join(
    ROOT, "code", "ecmodels-paper", "SupData_Scripts", "data", "ecYeast8_updated_lactate.xml"
)

# ---------------------------------------------------------------- parameters
PKA = 3.86          # lactic acid, Abbott et al. 2008
MW_GLUCOSE = 0.18016  # g mmol^-1
MW_LACTATE = 0.09008  # g mmol^-1
NGAM = 0.7          # mmol ATP gDW^-1 h^-1, ecYeast8 r_4046

# toxicity module (fitted, see calibrate.py)
TOX = {"k_h": None, "pK_phi": None, "k_a": None}

GENE = {
    "PDC1": "YLR044C", "PDC5": "YLR134W", "PDC6": "YGR087C",
    "GPD1": "YDL022W", "GPD2": "YOL059W",
    "ADH1": "YOL086C", "ADH2": "YMR303C", "ADH3": "YMR083W",
    "ADH4": "YGL256W", "ADH5": "YBR145W",
    "DLD1": "YDL174C", "JEN1": "YKL217W", "CYB2": "YML054C",
}


def load_model():
    m = cobra.io.read_sbml_model(MODEL_PATH)
    m.solver.configuration.verbosity = 0
    m.reactions.r_2111.bounds = (0.0, 1000.0)      # biomass cap in the distributed file
    m.reactions.r_1992_REV.bounds = (0.0, 1000.0)  # oxygen
    m.reactions.r_1714_REV.bounds = (0.0, 1000.0)  # glucose
    return m


def apply_genotype(model, deletions=(), ldh_copies=1.0):
    """Delete genes by ORF id and scale LDH capacity by copy number."""
    orfs = [GENE[g] for g in deletions if g in GENE]
    present = [g for g in orfs if g in model.genes]
    if present:
        cobra.manipulation.delete_model_genes(model, present)
    if ldh_copies != 1.0:
        r = model.reactions.draw_prot_ldh
        # higher copy number -> lower protein cost per unit flux
        for met, coeff in list(r.metabolites.items()):
            if met.id == "prot_pool":
                r.add_metabolites({met: coeff * (1.0 / ldh_copies) - coeff})
    return model


def speciation(total_lactate_mM, pH):
    f_ha = 1.0 / (1.0 + 10.0 ** (pH - PKA))
    return f_ha * total_lactate_mM, (1.0 - f_ha) * total_lactate_mM


def phi(pH, pK_phi):
    """Relative membrane permeability of the undissociated acid.

    pK_phi=None corresponds to a pH-independent permeability (the two-parameter
    model actually used here).
    """
    if pK_phi is None:
        return 1.0
    return 1.0 / (1.0 + 10.0 ** (pK_phi - pH))


def acid_maintenance(total_lactate_mM, pH, tox=None):
    """Extra ATP demand caused by lactic acid (mmol ATP gDW^-1 h^-1)."""
    t = tox or TOX
    ha, an = speciation(total_lactate_mM, pH)
    return t["k_h"] * phi(pH, t["pK_phi"]) * ha + t["k_a"] * an, ha, an


# ---------------------------------------------------------------- steady state
def chemostat_yield(model, D=0.1, maintenance=0.0, substrate="r_1714_REV"):
    """Biomass yield (g g^-1) in a glucose-limited chemostat at dilution rate D."""
    m = model.copy()
    m.reactions.r_2111.bounds = (D, D)
    m.reactions.r_4046.bounds = (NGAM + maintenance, NGAM + maintenance)
    m.objective = substrate
    m.objective.direction = "min"
    s = m.optimize()
    if s.status != "optimal":
        return np.nan
    q = s.fluxes[substrate]
    return D / (q * MW_GLUCOSE) if q > 0 else np.nan


def _fba_step(model, glucose_ub, maintenance):
    """One quasi-steady-state LP: maximise growth under a substrate cap."""
    m = model
    m.reactions.r_1714_REV.bounds = (0.0, max(glucose_ub, 1e-9))
    m.reactions.r_4046.bounds = (NGAM + maintenance, NGAM + maintenance)
    m.objective = "r_2111"
    s = m.optimize()
    if s.status != "optimal" or s.objective_value is None:
        return None
    return s


# ---------------------------------------------------------------- integrators
def simulate_batch(model, pH=6.0, X0=0.1, S0=44.4, V0=1.0, t_end=120.0, dt=0.1,
                   qs_max=10.0, Ks=0.28, tox=None, drift_buffer=None):
    """Batch culture. drift_buffer (mM per pH unit) enables uncontrolled pH."""
    m = model.copy()
    X, S, L, V, t = X0, S0, 0.0, V0, 0.0
    rows = []
    pH_now = pH
    while t <= t_end + 1e-9:
        if drift_buffer is not None:
            pH_now = max(2.2, pH - L / drift_buffer)
        mad, ha, an = acid_maintenance(L, pH_now, tox)
        q = qs_max * S / (Ks + S) if S > 0 else 0.0
        sol = _fba_step(m, q, mad)
        if sol is None:
            mu, ql, qg = 0.0, 0.0, 0.0
        else:
            mu = sol.objective_value
            ql = sol.fluxes["r_1546"]
            qg = sol.fluxes["r_1714_REV"]
        rows.append(dict(t=t, X=X, S=S, L=L, pH=pH_now, mu=mu, q_lac=ql,
                         q_glc=qg, m_acid=mad))
        X = max(0.0, X + mu * X * dt)
        S = max(0.0, S - qg * X * dt)
        L = max(0.0, L + ql * X * dt)
        t += dt
        if X <= 1e-6 and S <= 1e-6:
            break
    return pd.DataFrame(rows)


def simulate_fedbatch(model, pH=6.0, X0=0.5, S0=20.0, V0=0.5, t_end=100.0, dt=0.1,
                      qs_max=10.0, Ks=0.28, tox=None, S_feed=500.0,
                      mu_set=0.12, S_target=5.0, drift_buffer=None):
    """Exponential fed-batch with a feedback controller on residual glucose."""
    m = model.copy()
    X, S, L, V, t = X0, S0, 0.0, V0, 0.0
    rows = []
    F = 0.0
    pH_now = pH
    while t <= t_end + 1e-9:
        if drift_buffer is not None:
            pH_now = max(2.2, pH - L / drift_buffer)
        mad, ha, an = acid_maintenance(L, pH_now, tox)
        q = qs_max * S / (Ks + S) if S > 0 else 0.0
        sol = _fba_step(m, q, mad)
        if sol is None:
            mu, ql, qg = 0.0, 0.0, 0.0
        else:
            mu, ql, qg = sol.objective_value, sol.fluxes["r_1546"], sol.fluxes["r_1714_REV"]
        rows.append(dict(t=t, X=X, S=S, L=L, V=V, F=F, pH=pH_now, mu=mu, q_lac=ql, q_glc=qg, m_acid=mad))
        # feedback: keep residual glucose near S_target by adjusting dilution
        F = max(0.0, (mu * X * V + (S_target - S) * 0.3 * V) / max(S_feed - S, 1e-6))
        F = min(F, 0.02 * V)
        X = max(0.0, X + (mu * X - F / V * X) * dt)
        S = max(0.0, S + (-qg * X + F / V * (S_feed - S)) * dt)
        L = max(0.0, L + (ql * X - F / V * L) * dt)
        V += F * dt
        t += dt
        if V > 5 * V0:
            break
    return pd.DataFrame(rows)


def summarise(df):
    last = df.iloc[-1]
    titer = last["L"] * MW_LACTATE           # g/L
    glucose_used = 0.0
    prod = titer / max(last["t"], 1e-9)
    return dict(titer_gL=titer, time_h=last["t"], productivity=prod,
                biomass=last["X"], pH=last.get("pH", np.nan))
