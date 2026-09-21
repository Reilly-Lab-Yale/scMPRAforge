#!/usr/bin/env python
"""Poisson is ruled out: observed dispersion vs its simulated null.

Draws the result of overdispersion.py, which is the step that motivates an
overdispersed count family at all. For each cell type the Pearson dispersion of a Poisson fit conditioning on
CRE identity, phi = X2/(n-p), is plotted against the envelope of phi obtained
by refitting data simulated from that same Poisson fit. Under Poisson phi is
1; the simulated null pins down how far it can stray by chance at these
sample sizes and mean levels, where the asymptotic chi-square reference is
not calibrated.

Cell types are drawn as one strip rather than one labelled row each. Which
cell type carries which phi is not the claim -- the claim is that every one
of them lands decades away from the null -- and naming them invited the
reader to compare identities the panel cannot speak to.

    python analyses/model_selection/plot_overdispersion.py
"""
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

# Editable <text> so the manuscript sync pipeline's transforms apply, matching
# the other figures.
plt.rcParams["svg.fonttype"] = "none"

BASE = pathlib.Path(__file__).resolve().parent
OUT = BASE / "output"
TSV = BASE / "overdispersion.tsv"

BLUE, ORANGE = "#0072b2", "#d55e00"     # Okabe-Ito (house palette)
INK, MUTED = "#1a1a1a", "#6b6b6b"


def main():
    d = pd.read_csv(TSV, sep="\t")
    assert d.dataset.nunique() == 1, f"expected one dataset, got {sorted(d.dataset.unique())}"
    assert d.cell_type.is_unique, "cell types are not unique; one row per cell type expected"
    assert d.phi.notna().all(), "missing phi"
    # The whole point of the figure: observed dispersion sits outside the null.
    assert (d.phi > d.phi_null_max).all(), (
        "some cell type's observed phi is within its simulated null -- "
        "the figure would be claiming something the data does not support")
    assert (d.aic_pois > d.aic_nb).all(), "Poisson beats NB by AIC somewhere"

    d = d.sort_values("phi")
    # Included at 0.70\textwidth (498.66pt), so the canvas is drawn 4.83in
    # wide and point sizes below mean what they say on the page.
    fig, ax = plt.subplots(figsize=(4.83, 1.35))

    # Deterministic beeswarm: phi clusters between 15 and 20, and simply
    # alternating rows still collides there. Each point takes the row nearest
    # the axis that no placed point already occupies within MIN_SEP (measured
    # in log10 phi, the axis the eye actually reads).
    MIN_SEP, ROWS = 0.035, (0.0, 0.21, -0.21, 0.42, -0.42)
    placed, y = [], []
    for xi in d.phi:
        for row in ROWS:
            if all(abs(np.log10(xi) - np.log10(xj)) >= MIN_SEP
                   for xj, r in placed if r == row):
                break
        placed.append((xi, row))
        y.append(row)
    y = np.array(y)

    # Null envelope: the full spread of phi across all simulations, pooled
    # over cell types. It is narrow enough that per-cell-type bands would
    # overplot into a single stripe anyway.
    lo = float((d.phi_null_mean - 3 * d.phi_null_sd).min())
    hi = float(d.phi_null_max.max())
    ax.axvspan(lo, hi, color=MUTED, alpha=0.18, zorder=1)
    ax.axvline(1.0, color=MUTED, lw=1.0, zorder=2)

    # 2px surface ring so dots stay separable where they overlap.
    ax.scatter(d.phi, y, s=42, color=BLUE, zorder=4,
               edgecolor="white", linewidth=1.0)

    ax.set_yticks([])
    ax.set_ylim(-0.95, 0.55)
    ax.set_xscale("log")
    ax.set_xlim(0.6, float(d.phi.max()) * 1.25)
    ax.set_xlabel("Pearson dispersion $\\phi$ of the Poisson fit "
                  "(1 = Poisson)", fontsize=9, color=INK)

    med = float(d.phi.median())
    # The median is a summary, not one of the dots, so it gets a rule rather
    # than a leader pointing into the gap between two points.
    ax.vlines(med, -0.55, 0.45, color=MUTED, lw=0.9, ls=(0, (3, 2)), zorder=2)
    ax.text(med, -0.66, f"median {med:.1f}", ha="center", va="top",
            fontsize=7.5, color=INK)

    ax.xaxis.set_major_formatter(
        matplotlib.ticker.FuncFormatter(lambda v, _: f"{v:g}"))
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.grid(True, axis="x", color="#e6e6e6", lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    ax.spines["bottom"].set_color("#cccccc")
    ax.tick_params(colors=MUTED, labelsize=8, length=0)
    ax.annotate("Poisson null",
                xy=(hi, 0.30), xytext=(1.35, 0.42),
                fontsize=7.5, color=MUTED, va="center", ha="left",
                arrowprops=dict(arrowstyle="-", color=MUTED, lw=0.8,
                                shrinkA=0, shrinkB=2))

    print(f"phi: {d.phi.min():.1f} to {d.phi.max():.1f} (median {d.phi.median():.1f})")
    print(f"null: mean {d.phi_null_mean.mean():.2f}, max {hi:.2f}")

    fig.tight_layout()
    OUT.mkdir(exist_ok=True)
    for ext in ("svg", "png"):
        fig.savefig(OUT / f"overdispersion_shendure.{ext}", dpi=200, bbox_inches="tight")
    plt.close(fig)
    print(f"wrote {OUT/'overdispersion_shendure.svg'} and .png")


if __name__ == "__main__":
    main()
