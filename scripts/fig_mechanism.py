"""Conceptual figure: the two-mechanism model of lactic-acid toxicity."""

import os

import matplotlib
matplotlib.use("Agg")
from PIL import Image
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch, Rectangle

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIG = os.path.join(ROOT, "figures")

C_TEXT = "#1a1a1a"
C_MED = "#FDF3E3"
C_CELL = "#EAF2FB"
C_MEM = "#4D4D4D"
C_HA = "#A33B1C"
C_AN = "#0F6E56"
C_ENZ = "#7A4A12"
C_ATP = "#4B439B"

plt.rcParams.update({"font.family": "serif",
                     "font.serif": ["Times New Roman"], "font.size": 7,
                     "mathtext.fontset": "custom",
                     "mathtext.rm": "Times New Roman",
                     "mathtext.it": "Times New Roman:italic",
                     "mathtext.bf": "Times New Roman:bold",
                     "svg.fonttype": "none"})


def arrow(ax, p1, p2, color, lw=1.1, ls="-", rad=0.0, style="-|>"):
    ax.add_patch(FancyArrowPatch(p1, p2, arrowstyle=style, mutation_scale=7,
                                 color=color, lw=lw, linestyle=ls,
                                 connectionstyle=f"arc3,rad={rad}", shrinkA=1.5, shrinkB=1.5))


def box(ax, x, y, w, h, label, edge, fs=7, fc="white"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.25,rounding_size=1.0",
                                facecolor=fc, edgecolor=edge, lw=0.8, zorder=4))
    ax.text(x + w / 2, y + h / 2, label, ha="center", va="center", fontsize=fs, color=C_TEXT, zorder=5)



def _flatten(path, dpi=900):
    """Strip the alpha channel so the TIFF is RGB, as print workflows expect."""
    im = Image.open(path)
    if im.mode != "RGB":
        rgb = Image.new("RGB", im.size, "white")
        rgb.paste(im, mask=im.split()[-1] if im.mode in ("RGBA", "LA") else None)
        rgb.save(path, dpi=(dpi, dpi), compression="tiff_lzw")

def main():
    fig, ax = plt.subplots(figsize=(5.9, 3.4))
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 58)
    ax.axis("off")

    ax.add_patch(Rectangle((0, 36.5), 100, 21.5, facecolor=C_MED, edgecolor="none", zorder=0))
    ax.add_patch(FancyBboxPatch((2, 2), 96, 34, boxstyle="round,pad=0.3,rounding_size=2.5",
                                facecolor=C_CELL, edgecolor="#B5D4F4", lw=0.6, zorder=1))

    for y in (34.4, 35.4):
        ax.plot([2, 98], [y, y], color=C_MEM, lw=1.3, zorder=3, solid_capstyle="butt")

    ax.text(3, 55.4, "Extracellular, pH 4.5", fontsize=7, color="#5F5E5A")
    ax.text(4, 8.2, "Cytosol, pH 6.8", fontsize=7, color="#185FA5")

    # ---------------- outside: the two species
    ax.text(7.5, 46.6, "undissociated acid (HA)", fontsize=7, color=C_HA)
    for x in (13, 17.5, 22):
        ax.add_patch(Circle((x, 43.2), 1.0, facecolor=C_HA, edgecolor="none", alpha=0.85, zorder=5))
    arrow(ax, (17.5, 41.9), (17.5, 36.0), C_HA, lw=1.3)
    ax.text(18.6, 38.2, "passive diffusion", fontsize=6.4, color=C_HA)

    ax.text(38, 46.6, "lactate anion (A$^-$)", fontsize=7, color=C_AN)
    for x in (40, 44.5):
        ax.add_patch(Circle((x, 43.2), 1.0, facecolor=C_AN, edgecolor="none", alpha=0.85, zorder=5))
    arrow(ax, (42.2, 41.9), (42.2, 37.2), C_AN, lw=1.0, ls=(0, (2, 1.5)))
    ax.plot([40.9, 43.6], [38.6, 35.9], color=C_AN, lw=1.1, zorder=6)
    ax.plot([40.9, 43.6], [35.9, 38.6], color=C_AN, lw=1.1, zorder=6)

    # ---------------- inside: dissociation
    ax.text(24, 30.6, "HA  \u2192  H$^+$  +  A$^-$", fontsize=8, color=C_HA, ha="center")
    ax.text(24, 27.6, "acidification of the cytosol", fontsize=6.2, color="#8C3B1E", ha="center")

    # ---------------- PMA1 in the membrane
    ax.add_patch(Rectangle((52, 33.4), 9, 3.0, facecolor="#F7E4E4", edgecolor=C_ENZ, lw=0.9, zorder=6))
    ax.text(56.5, 34.9, "PMA1", ha="center", va="center", fontsize=6.4, color=C_ENZ, zorder=7)
    arrow(ax, (56.5, 31.4), (56.5, 33.2), C_ATP, lw=1.2)
    arrow(ax, (56.5, 36.6), (56.5, 39.4), C_ATP, lw=1.2)
    ax.text(57.6, 37.6, "H$^+$", fontsize=6.8, color=C_ATP)
    ax.text(55.6, 30.0, "ATP", fontsize=6.4, color=C_ATP, ha="right")

    ax.text(34.5, 23.4, "proton export\n(k$_H$)", fontsize=6.4, color=C_ATP, ha="center")
    arrow(ax, (36.5, 25.6), (51.0, 31.6), C_ATP, lw=0.9, ls=(0, (1.5, 1.5)), rad=-0.20)

    # ---------------- central carbon pathway
    box(ax, 6, 12, 15, 5.6, "glucose", "#888780")
    box(ax, 27, 12, 13, 5.6, "pyruvate", "#888780")
    box(ax, 45, 12, 15, 5.6, "D-lactate", C_AN)
    arrow(ax, (21.3, 14.8), (26.7, 14.8), "#5F5E5A")
    arrow(ax, (40.3, 14.8), (44.7, 14.8), "#5F5E5A")
    ax.text(42.5, 16.4, "LDH", fontsize=6.4, color=C_ENZ, ha="center")
    arrow(ax, (12.5, 36.0), (12.5, 17.8), "#888780", lw=1.0)
    arrow(ax, (60.3, 14.8), (65.0, 14.8), C_AN, lw=1.0)
    ax.text(66.0, 14.8, "export", fontsize=6.4, color=C_AN, va="center")

    # ---------------- anion fate
    ax.text(35.5, 6.4, "A$^-$ retained: osmotic load and iron chelation (k$_A$)",
            fontsize=6.4, color=C_AN, ha="center")

    # ---------------- ATP ledger
    ax.add_patch(FancyBboxPatch((67, 2.4), 30, 10.4, boxstyle="round,pad=0.3,rounding_size=1.2",
                                facecolor="white", edgecolor="#B4B2A9", lw=0.6, zorder=4))
    ax.text(68.6, 11.3, "both terms enter as", fontsize=6.2, color="#5F5E5A", zorder=5)
    ax.text(68.6, 8.8, "m$_{acid}$ = k$_H$[HA] + k$_A$[A$^-$]", fontsize=7.0, color=C_TEXT, zorder=5)
    ax.text(68.6, 5.9, "an ATP demand in the", fontsize=6.2, color="#5F5E5A", zorder=5)
    ax.text(68.6, 4.1, "enzyme-constrained model", fontsize=6.2, color="#5F5E5A", zorder=5)

    # ---------------- consequence
    ax.text(79.5, 20.8, "ATP diverted from growth\nand product formation",
            fontsize=6.5, color=C_ATP, ha="left")
    arrow(ax, (61.0, 29.4), (78.5, 23.4), C_ATP, lw=0.9, ls=(0, (1.5, 1.5)), rad=-0.18)

    fig.tight_layout()
    stem = os.path.join(FIG, "fig1_mechanism")
    fig.savefig(stem + ".png", bbox_inches="tight")
    fig.savefig(stem + ".svg", bbox_inches="tight")
    fig.savefig(stem + ".tif", dpi=900, bbox_inches="tight",
                pil_kwargs={"compression": "tiff_lzw"}, facecolor="white")
    _flatten(stem + ".tif")
    plt.close(fig)
    print("saved fig1_mechanism .png/.svg/.tif")


if __name__ == "__main__":
    main()
