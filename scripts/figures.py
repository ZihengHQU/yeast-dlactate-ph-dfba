"""Publication figures for the manuscript."""

import json
import os

import matplotlib
matplotlib.use("Agg")
from PIL import Image
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RES = os.path.join(ROOT, "results")
FIG = os.path.join(ROOT, "figures")
os.makedirs(FIG, exist_ok=True)

plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["Times New Roman"],
    "mathtext.fontset": "custom",
    "mathtext.rm": "Times New Roman",
    "mathtext.it": "Times New Roman:italic",
    "mathtext.bf": "Times New Roman:bold",
    "font.size": 8,
    "axes.linewidth": 0.6,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "savefig.dpi": 400,
    "svg.fonttype": "none",
})

C = {"ec": "#b03a2e", "y8": "#2c3e50", "fit": "#1f618d", "lit": "#7f8c8d",
     "alt": "#117a65", "warn": "#b7950b"}


def save(fig, name):
    """Write PNG (preview), SVG (vector) and 900 dpi TIFF (print)."""
    stem = os.path.splitext(name)[0]
    fig.savefig(os.path.join(FIG, stem + ".png"), bbox_inches="tight")
    fig.savefig(os.path.join(FIG, stem + ".svg"), bbox_inches="tight")
    fig.savefig(os.path.join(FIG, stem + ".tif"), dpi=900, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"}, facecolor="white")
    _flatten(os.path.join(FIG, stem + ".tif"))
    plt.close(fig)
    print("saved", stem, ".png/.svg/.tif")



def _flatten(path, dpi=900):
    """Strip the alpha channel so the TIFF is RGB, as print workflows expect."""
    im = Image.open(path)
    if im.mode != "RGB":
        rgb = Image.new("RGB", im.size, "white")
        rgb.paste(im, mask=im.split()[-1] if im.mode in ("RGBA", "LA") else None)
        rgb.save(path, dpi=(dpi, dpi), compression="tiff_lzw")

def fig1():
    p = json.load(open(os.path.join(ROOT, "data", "toxicity_params.json")))
    m50, k_h, k_a = p["m50"], p["k_h"], p["k_a"]
    pka = 3.86
    obs = [(3.0, 500), (3.5, 750), (5.0, 900)]
    phs = np.linspace(2.6, 5.6, 200)

    def c_required(pH, kk_h, kk_a):
        f = 1 / (1 + 10 ** (pH - pka))
        return m50 / (kk_h * f + kk_a * (1 - f))

    k_pure = m50 / (900 * (1 / (1 + 10 ** (5.0 - pka))))

    fig, ax = plt.subplots(1, 2, figsize=(6.7, 2.75))
    ax[0].plot(phs, [c_required(x, k_pure, 0) for x in phs], "--", color=C["lit"],
               lw=1.2, label="undissociated acid only")
    ax[0].plot(phs, [c_required(x, k_h, k_a) for x in phs], "-", color=C["ec"],
               lw=1.6, label="two-mechanism model")
    ax[0].scatter([o[0] for o in obs], [o[1] for o in obs], s=32, color="k", zorder=5,
                  label="Abbott 2008 chemostat")
    ax[0].set_yscale("log")
    ax[0].set_xlabel("Culture pH")
    ax[0].set_ylabel("Total lactic acid for 50% yield loss (mM)")
    ax[0].legend(frameon=False, fontsize=6.5, loc="upper left", handlelength=1.6, borderaxespad=0.2, bbox_to_anchor=(0.03, 0.94))
    ax[0].annotate("8.6x", xy=(3.46, 95), xytext=(2.95, 150), fontsize=7, color=C["lit"],
                   arrowprops=dict(arrowstyle="-", lw=0.5, color=C["lit"]))
    ax[0].annotate("1.35x", xy=(3.52, 575), xytext=(3.92, 800), fontsize=7, color=C["ec"],
                   arrowprops=dict(arrowstyle="-", lw=0.5, color=C["ec"]))
    ax[0].set_title("A", loc="left", fontsize=9, fontweight="bold")

    # panel B: predicted maintenance demand at the observed concentrations
    xs, ys = [], []
    for pH, conc in obs:
        f = 1 / (1 + 10 ** (pH - pka))
        xs.append(pH)
        ys.append((k_h * f * conc + k_a * (1 - f) * conc) / m50)
    ax[1].bar(range(3), ys, color=[C["ec"], C["alt"], C["fit"]], width=0.55)
    ax[1].axhline(1.0, color="k", lw=0.8, ls=":")
    ax[1].set_xticks(range(3))
    ax[1].set_xticklabels([f"pH {p}" for p, _ in obs])
    ax[1].set_ylabel("Predicted maintenance / m$_{50}$")
    ax[1].set_ylim(0, 1.45)
    ax[1].set_title("B", loc="left", fontsize=9, fontweight="bold")
    for a in ax:
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save(fig, "fig1_toxicity_model.png")


def fig2():
    f1 = os.path.join(RES, "traj_Baek_batch.csv")
    f2 = os.path.join(RES, "traj_Ishida_batch.csv")
    if not (os.path.exists(f1) and os.path.exists(f2)):
        print("fig2 skipped")
        return
    a = pd.read_csv(f1)
    b = pd.read_csv(f2)
    fig, ax = plt.subplots(1, 3, figsize=(6.7, 2.45))
    ax[0].plot(a["t"], a["L"] * 0.09008, color=C["ec"], lw=1.4, label="pH uncontrolled")
    ax[0].plot(b["t"], b["L"] * 0.09008, color=C["fit"], lw=1.4, label="pH 6.5")
    ax[0].set_xlabel("Time (h)")
    ax[0].set_ylabel("D-lactate (g L$^{-1}$)")
    ax[0].legend(frameon=False, fontsize=7)
    ax[1].plot(a["t"], a["S"] * 0.18016, color=C["ec"], lw=1.2)
    ax[1].plot(b["t"], b["S"] * 0.18016, color=C["fit"], lw=1.2)
    ax[1].set_xlabel("Time (h)")
    ax[1].set_ylabel("Glucose (g L$^{-1}$)")
    ax[2].plot(a["t"], a["pH"], color=C["ec"], lw=1.2)
    ax[2].set_xlabel("Time (h)")
    ax[2].set_ylabel("Culture pH")
    ax[2].set_ylim(2.0, 7.0)
    for i, x in enumerate(ax):
        x.spines[["top", "right"]].set_visible(False)
        x.set_title("ABC"[i], loc="left", fontsize=9, fontweight="bold")
    fig.tight_layout()
    save(fig, "fig2_batch_dynamics.png")


def fig3():
    f = os.path.join(RES, "counterfactual.csv")
    if not os.path.exists(f):
        print("fig3 skipped")
        return
    d = pd.read_csv(f)
    fig, ax = plt.subplots(figsize=(4.6, 4.1))
    lim = [5, 1200]
    ax.plot(lim, lim, "--", color=C["lit"], lw=0.8, zorder=1)
    ax.fill_between(lim, [x * 0.8 for x in lim], [x * 1.2 for x in lim],
                    color=C["lit"], alpha=0.12, lw=0, zorder=0)
    two = d[d.model == "two"]
    one = d[d.model == "single"]
    ax.scatter(two.titer_lit, two.titer, s=42, color=C["ec"], zorder=5,
               edgecolor="white", linewidth=0.4, label="two mechanisms")
    ax.scatter(one.titer_lit, one.titer, s=36, facecolor="white", edgecolor=C["fit"],
               linewidth=1.0, marker="s", zorder=4, label="undissociated acid only")
    ax.annotate("no predicted ceiling\nat pH 6.0", (112, 835), textcoords="offset points",
                xytext=(9, -16), fontsize=7, color=C["fit"])
    ax.annotate("9-fold low\nat pH 3.5", (82.6, 8.8), textcoords="offset points",
                xytext=(8, -2), fontsize=7, color=C["fit"])
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlim(lim); ax.set_ylim(lim)
    ax.set_xlabel("Reported titre (g L$^{-1}$)")
    ax.set_ylabel("Predicted titre (g L$^{-1}$)")
    ax.legend(frameon=False, fontsize=7.5, loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save(fig, "fig3_validation.png")


def fig4():
    f = os.path.join(RES, "ph_sweep.csv")
    g = os.path.join(RES, "genotypes.csv")
    if not os.path.exists(f):
        print("fig4 skipped")
        return
    d = pd.read_csv(f)
    fig, ax = plt.subplots(1, 2, figsize=(6.7, 2.95),
                           gridspec_kw={"width_ratios": [1.0, 1.15], "wspace": 0.55})
    ax[0].plot(d["pH"], d["titer"], "o-", color=C["ec"], lw=1.5, ms=3.5)
    ax[0].axvline(3.86, color=C["lit"], ls=":", lw=0.8)
    ax[0].text(3.92, d["titer"].min() + 1.5, "pKa", fontsize=6.5, color=C["lit"])
    ax[0].set_xlabel("Controlled pH")
    ax[0].set_ylabel("Final D-lactate (g L$^{-1}$)")
    ax[0].set_ylim(d["titer"].min() - 3, d["titer"].max() + 3)
    ax[0].set_title("A", loc="left", fontsize=9, fontweight="bold")

    if os.path.exists(g):
        e = pd.read_csv(g)
        e["short"] = e["genotype"].str.replace("pdc1Δ + gpdΔ + dld1Δ + jen1Δ", "+gpdΔ dld1Δ jen1Δ") \
                                 .str.replace("pdc1Δ + gpd1/2Δ", "+gpd1/2Δ")
        piv = e.pivot(index="short", columns="pH", values="titer").loc[
            ["pdc1Δ", "+gpd1/2Δ", "+gpdΔ dld1Δ jen1Δ", "full (Baek 2016)"]]
        piv.plot(kind="barh", ax=ax[1], color=[C["fit"], C["ec"]], width=0.72, legend=True)
        ax[1].set_xlabel("Final D-lactate (g L$^{-1}$)")
        ax[1].set_ylabel("")
        ax[1].tick_params(axis="y", labelsize=7)
        ax[1].set_xlim(0, 100)
        ax[1].legend(frameon=False, fontsize=6.5, title="pH", title_fontsize=6.5,
                     loc="center right")

        ax[1].set_title("B", loc="left", fontsize=9, fontweight="bold")
    for a in ax:
        a.spines[["top", "right"]].set_visible(False)
    fig.subplots_adjust(bottom=0.2, top=0.9)
    save(fig, "fig4_design_space.png")


def fig5():
    f = os.path.join(RES, "sensitivity.csv")
    if not os.path.exists(f):
        print("fig5 skipped")
        return
    d = pd.read_csv(f)
    fig, ax = plt.subplots(1, 2, figsize=(6.7, 2.7), sharey=False)
    for i, p in enumerate(["q_lac_max", "q_O2_max"]):
        sub = d[d["param"] == p]
        for pH, col in [(6.0, C["fit"]), (3.5, C["ec"])]:
            s2 = sub[sub["pH"] == pH].sort_values("value")
            ax[i].plot(s2["value"], s2["titer"], "o-", ms=3, lw=1.2, color=col, label=f"pH {pH}")
        ax[i].set_xlabel("q$_{lac,max}$ (mmol gDW$^{-1}$ h$^{-1}$)" if p == "q_lac_max"
                         else "q$_{O2,max}$ (mmol gDW$^{-1}$ h$^{-1}$)")
        ax[i].set_ylabel("Final D-lactate (g L$^{-1}$)")
        ax[i].legend(frameon=False, fontsize=7, loc="lower right")
        vals = sorted(sub["value"].unique())
        ax[i].set_xticks(vals)
        ax[i].set_xticklabels([f"{v:g}" for v in vals])
        ax[i].set_title("AB"[i], loc="left", fontsize=9, fontweight="bold")
        ax[i].spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    save(fig, "fig5_sensitivity.png")


if __name__ == "__main__":
    fig1()
    fig2()
    fig3()
    fig4()
    fig5()
