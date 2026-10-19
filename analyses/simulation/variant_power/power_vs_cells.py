#!/usr/bin/env python3
"""Pairwise power at the published median variant effects, against cell number.

Manuscript Fig 3C: one panel per dataset, under that dataset's Fig 3B heatmap.
Each point is one cell type. Its y is the power of the heatmap for that cell
type read along the two Kircher et al. reference lines (|log2 FC| 0.263 and
0.398, interpolated between grid rows) and averaged over the baseline columns,
which are log-spaced from the 5th to the 95th percentile of the dataset's
observed activity. Its x is the number of cells of that type.

Datasets are kept in separate panels because cell number is not the only thing
that differs between them (MOI, library, baseline activity, dispersion all do).
Within one dataset the simulations share all of those, and each cell type
carries the cell count and CRE set observed in it.

    python power_vs_cells.py      # on Bouchet; needs scMPRAforge for the cell counts

Writes output/power_vs_cells/power_vs_cells_<dataset>.svg and
output/power_vs_cells/power_vs_cells.tsv.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
sys.path.insert(0, str(HERE))

import scMPRAforge as scm
import kircher_reference as kr
import panel_style as ps

OUT = HERE / "output" / "power_vs_cells"

BOUNDS = {
    "shendure": scm.SHENDURE_BOUNDS,
    "cohen": scm.COHEN_BOUNDS,
    "seelig": scm.SEELIG_BOUNDS,
}

# The two series wear the same colour and marker as their reference lines on
# the Fig 3B heatmaps, so the eye carries a line down from B to its series in C.
SERIES = [
    ("activating", kr.ACTIVATING_LOG2FC, kr._ACTIVATING_STYLE, "activating"),
    ("repressing", kr.REPRESSING_LOG2FC, kr._REPRESSING_STYLE, "repressing"),
]

# One x range for all three panels, so a cell count sits at the same place in
# each. Spans the smallest (680, Lalanne cardiomyocytes) to the largest
# (18,633, Zhao rods) cell type with a margin either side.
X_LIM = (450, 30000)
X_TICKS = [1000, 3000, 10000]

# Room for the dataset title only.
TOP_IN = 0.24
FIG_H = TOP_IN + ps.PLOT_H_IN + ps.BOTTOM_IN


def summarise(dataset):
    """One row per cell type: cells, and mean power at each reference effect."""
    df = pd.read_parquet(HERE / dataset / "output" / f"{dataset}_pairwise_power_df.parquet")
    grid = (df.groupby(["cell_type", "baseline_mu", "log2_fc"])["reject_null"]
              .mean().rename("power").reset_index())
    cells = BOUNDS[dataset].cells_per_cell_type
    assert set(grid["cell_type"]) == set(cells.index), (
        f"{dataset}: cell types differ between power grid "
        f"{sorted(grid['cell_type'].unique())} and bounds {sorted(cells.index)}")
    assert len(cells) == ps.N_CELL_TYPES[dataset], (
        f"{dataset}: {len(cells)} cell types, expected {ps.N_CELL_TYPES[dataset]}")

    rows = []
    for ct, g in grid.groupby("cell_type"):
        pivot = g.pivot(index="log2_fc", columns="baseline_mu", values="power").sort_index()
        assert not pivot.isna().any().any(), f"{dataset}/{ct}: holes in the power grid"
        fc = pivot.index.to_numpy()
        row = {"dataset": dataset, "cell_type": ct,
               "cell_type_display": ps.display_ct(ct, dataset),
               "cells": int(cells[ct]), "n_baselines": pivot.shape[1]}
        for name, target, _, _ in SERIES:
            assert fc.min() <= target <= fc.max(), (
                f"{dataset}/{ct}: {target} outside grid {fc.min()}-{fc.max()}")
            per_baseline = np.array([np.interp(target, fc, pivot[c]) for c in pivot.columns])
            row[f"power_{name}"] = per_baseline.mean()
            row[f"power_{name}_max"] = per_baseline.max()
        rows.append(row)
    out = pd.DataFrame(rows).sort_values("cells").reset_index(drop=True)
    assert out["cells"].gt(0).all(), f"{dataset}: a cell type with no cells"
    # Domain check: the larger effect is never harder to detect.
    assert (out["power_repressing"] >= out["power_activating"] - 0.02).all(), (
        f"{dataset}: power at |log2 FC| 0.40 below power at 0.26:\n{out}")
    return out


def plot(summary, dataset, legend):
    fig, ax = plt.subplots(figsize=(ps.FIG_W, FIG_H))

    ax.axhline(0.8, color=ps.MUTED, lw=0.6, ls=(0, (3, 2)), zorder=1)
    for name, _, style, label in SERIES:
        ax.plot(summary["cells"], summary[f"power_{name}"],
                color=style["color"], lw=0.8, zorder=2)
        ax.plot(summary["cells"], summary[f"power_{name}"],
                linestyle="none", marker=style["marker"], ms=3.6,
                color=style["color"], markeredgecolor="white",
                markeredgewidth=0.5, zorder=3, label=label)

    # Ring the cell type Fig 3B shows. The caption says which it is: a label in
    # the plot collides with the 80% rule or the neighbouring points.
    ref = summary[summary["cell_type"] == "reference"]
    assert len(ref) == 1, f"{dataset}: {len(ref)} reference rows"
    ref = ref.iloc[0]
    for name, _, _, _ in SERIES:
        ax.plot(ref["cells"], ref[f"power_{name}"], linestyle="none", marker="o",
                ms=7.5, markerfacecolor="none", markeredgecolor=ps.INK,
                markeredgewidth=0.7, zorder=4)

    ax.set_xscale("log")
    ax.set_xlim(*X_LIM)
    ax.set_xticks(X_TICKS)
    ax.set_xticklabels([f"{t:,}" for t in X_TICKS])
    ax.xaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.set_ylim(0, 1.04)
    ax.set_yticks([0, 0.5, 1])
    ax.set_yticklabels(["0", "0.5", "1"])
    ax.set_xlabel("cells", color=ps.INK, labelpad=2)
    ax.set_ylabel("power", color=ps.INK, labelpad=2)
    ax.tick_params(colors=ps.INK, pad=2, length=2, width=0.6, which="both")
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_linewidth(0.6)

    if legend:
        # Repressing first: its series sits above activating in every panel,
        # so the key reads in the same top-to-bottom order as the plot.
        handles, labels = ax.get_legend_handles_labels()
        assert labels == [lab for *_, lab in SERIES], labels
        ax.legend(handles[::-1], labels[::-1], loc="lower right", frameon=False, handlelength=1.6,
                  handletextpad=0.4, borderaxespad=0.1, labelspacing=0.3)

    x = ps.LEFT_IN / ps.FIG_W
    texts = [fig.text(x, 1 - ps.TITLE_TOP_IN / FIG_H, ps.DISPLAY_NAME[dataset],
                      ha="left", va="center", fontsize=9, color=ps.INK)]
    fig.subplots_adjust(left=ps.LEFT_IN / ps.FIG_W, right=1 - ps.RIGHT_IN / ps.FIG_W,
                        top=1 - TOP_IN / FIG_H, bottom=ps.BOTTOM_IN / FIG_H)
    ps._assert_fits(fig, texts)
    path = OUT / f"power_vs_cells_{dataset}.svg"
    # Fixed canvas, no tight bbox, so the plot areas line up with Fig 3B above.
    fig.savefig(path, format="svg")
    plt.close(fig)
    print(f"  wrote {path}")


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    tables = []
    for i, dataset in enumerate(BOUNDS):
        s = summarise(dataset)
        assert s["cells"].between(*X_LIM).all(), (
            f"{dataset}: cell counts {s['cells'].min()}-{s['cells'].max()} outside {X_LIM}")
        print(f"{ps.DISPLAY_NAME[dataset]}")
        print(s[["cell_type_display", "cells", "power_activating", "power_repressing"]]
              .round(2).to_string(index=False))
        plot(s, dataset, legend=(i == 0))
        tables.append(s)
    pd.concat(tables).to_csv(OUT / "power_vs_cells.tsv", sep="\t", index=False,
                             float_format="%.4f")


if __name__ == "__main__":
    main()
