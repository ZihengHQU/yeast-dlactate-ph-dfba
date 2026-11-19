# pH-resolved enzyme-constrained dFBA for D-lactate production in *Saccharomyces cerevisiae*

Code and data for the manuscript *"Two mechanisms of lactic-acid toxicity explain pH-dependent
D-lactate tolerance in Saccharomyces cerevisiae: a pH-resolved enzyme-constrained dynamic flux
balance analysis"*.

The framework replaces the standard single-species description of weak-acid toxicity with two
additive energetic penalties — one for the undissociated acid, one for the anion — and embeds
them in enzyme-constrained dynamic flux balance analysis. The two coefficients are calibrated on
published chemostat physiology and then applied without further fitting.

## Key numbers

| | two mechanisms | undissociated acid only |
|---|---|---|
| mean absolute error, 8 strains | 33.4% | 116% |
| Spearman ρ | 0.51 | 0.03 |
| predicted titre, pH 3.5 | 55.4 g L⁻¹ | 8.8 g L⁻¹ (observed 82.6) |
| predicted titre, pH 6.0 | 95.6 g L⁻¹ | 835 g L⁻¹ (no ceiling) |

## Repository layout

```
src/core.py              model construction, toxicity module, dFBA integrators
scripts/calibrate.py     fits k_H and k_A on the Abbott 2008 chemostat data
scripts/run_all.py       pH sweep, genotype comparison, capacity sensitivity
scripts/counterfactual.py  eight-strain validation with the single-mechanism control
scripts/param_uncertainty.py  propagation of ±20% in the two coefficients
scripts/figures.py       data figures (PNG / SVG / 900 dpi TIFF)
scripts/fig_mechanism.py conceptual Figure 1
scripts/make_source_data.py  per-figure source data + MANIFEST
data/                    parameters, literature dataset, calibrated coefficients, references
results/                 simulation output (CSV)
source_data/             per-figure source data
figures/                 figures (PNG and SVG; TIFF generated on demand)
```

## Environment

Python 3.13 with:

```
cobra==0.32.1
pandas==2.3.3
numpy==2.5.5
scipy==1.18.1
matplotlib==3.11.2
python-docx
pillow
pypdfium2
```

Gurobi 13.0 was used as the LP solver (a GLPK fallback works but is 20–50× slower).
`gurobipy` must match the installed licence version.

## Model file

The enzyme-constrained model is **not** redistributed here because of its size. Download
`ecYeast8_updated_lactate.xml` from the supplementary material of

> Moreno-Paz S, Schmitz J, Martins dos Santos VAP, Suarez-Diez M. Enzyme-constrained models
> predict the dynamics of *Saccharomyces cerevisiae* growth in continuous, batch and fed-batch
> bioreactors. Microb Biotechnol. 2022;15:1435-1445.

and place it at `code/ecmodels-paper/SupData_Scripts/data/ecYeast8_updated_lactate.xml`,
which is where `src/core.py` looks for it (`MODEL_PATH`).

**Note** that the distributed file caps the biomass reaction at 0.13 h⁻¹; `load_model()` raises
this bound. Forgetting it produces a spurious growth ceiling that is easily mistaken for a
protein-pool limit.

## Reproduce

```
python scripts/calibrate.py          # ~2 min  -> data/toxicity_params.json
python scripts/counterfactual.py     # ~20 min -> results/counterfactual.csv
python scripts/run_all.py            # ~45 min -> results/*.csv
python scripts/figures.py            # figures as PNG, SVG and 900 dpi TIFF
python scripts/make_source_data.py   # source_data/ + MANIFEST.csv
```

## Model summary

Toxicity enters as an additional ATP demand on the maintenance reaction of ecYeast8:

```
m_acid = k_H · [HA] + k_A · [A⁻]
```

with `k_H = 0.0453` and `k_A = 0.0220` mmol ATP per mmol, obtained from the two controlled-pH
chemostat conditions of Abbott et al. (2008). The pH 3.5 condition of the same study is held out
and used as a test.

## Data sources

Validation strains come from Ishida 2006, Baek 2016, Baek 2017, Zhong 2019, Rahmasari 2025 and
Sornlek 2022; see `data/literature_strains.csv` for genotypes, conditions and DOIs. The chemostat
inhibition data are from Abbott et al. (2008), *Appl Environ Microbiol* 74:5759-5768.

## Contact

Ziheng Lee ([@ZihengHQU](https://github.com/ZihengHQU))
