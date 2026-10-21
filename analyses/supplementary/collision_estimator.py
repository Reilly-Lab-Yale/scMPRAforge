#!/usr/bin/env python3
"""Validate the occupancy estimator of collisions against simulated ground truth.

Manuscript Fig S9A. Each simulation delivers MPRA barcodes into cells and
counts collisions two ways:

  ground truth  deliveries minus distinct (cell, barcode) pairs, which needs
                to know every delivery -- including the repeats sequencing
                cannot see
  estimated     the occupancy estimator (Methods, Barcode collisions), which
                sees only the distinct barcodes in each cell: a cell holding
                k distinct of B library barcodes is inferred to have received
                log(1 - k/B) / log(1 - 1/B) deliveries. B is itself estimated
                as the number of distinct barcodes observed anywhere.

The grid crosses total cells, MOI and library size (475 simulations). The
per-cell delivery count is drawn from the Lalanne et al. transfection model
with its mean set to the MOI, as in the simulation framework.

Cells and barcodes are integer IDs rather than random sequences: random
10-mers for up to 801,000 cells can repeat, which would merge two cells.

    python collision_estimator.py simulate   # writes collision_estimator.tsv
    python collision_estimator.py plot       # writes coupon_collector_regression.svg

simulate needs scMPRAforge (for the fitted transfection model) and about
8 GB at the largest grid point (~76M deliveries); plot needs only the TSV.
"""
import sys
from itertools import product
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
TSV = HERE / "collision_estimator.tsv"
SVG = HERE / "coupon_collector_regression.svg"

TOTAL_CELLS = np.arange(1_000, 1_000_000, 200_000)
MOI = np.arange(5, 100, 5)
LIBRARY_SIZE = np.arange(10_000, 100_000, 20_000)
SEED = 20260224


def simulate_one(transfection_model, total_cells, library_size, rng):
    """Ground-truth and estimated collisions for one simulated replicate."""
    per_cell = transfection_model.draw_nb(int(total_cells), rng=rng)
    cell = np.repeat(np.arange(total_cells, dtype=np.int64), per_cell)
    bc = rng.integers(0, library_size, size=cell.size, dtype=np.int64)
    n = cell.size

    pairs = np.unique(cell * library_size + bc)          # distinct (cell, barcode)
    gt = n - pairs.size

    k = np.bincount(pairs // library_size, minlength=int(total_cells))
    b_hat = np.unique(bc).size
    assert b_hat <= library_size, f"{b_hat} distinct barcodes from a library of {library_size}"
    assert (k < b_hat).all(), "a cell holds every observed barcode; the estimator diverges"
    received = np.log1p(-k / b_hat) / np.log1p(-1 / b_hat)
    est = received.sum() - k.sum()
    return n, gt, est


def simulate():
    import scMPRAforge as scm
    params = pd.DataFrame(product(TOTAL_CELLS, MOI, LIBRARY_SIZE),
                          columns=["total_cells", "moi", "library_size"])
    assert len(params) == 475, f"grid has {len(params)} points, expected 475"
    rows = []
    for i, r in params.iterrows():
        bounds = scm.SHENDURE_BOUNDS.copy()
        bounds.set_effective_moi(r.moi)
        rng = np.random.default_rng(SEED + i)
        n, gt, est = simulate_one(bounds.transfection_model, r.total_cells,
                                  r.library_size, rng)
        rows.append((n, gt, est))
        print(f"{i + 1:3d}/475 cells={r.total_cells} moi={r.moi} "
              f"lib={r.library_size} deliveries={n} gt={gt} est={est:.1f}",
              flush=True)
    params[["deliveries", "gt_collisions", "estimated_collisions"]] = rows
    assert (params.gt_collisions >= 0).all(), "negative ground-truth collisions"
    assert (params.gt_collisions <= params.deliveries).all(), "more collisions than deliveries"
    params.to_csv(TSV, sep="\t", index=False)
    print(f"wrote {TSV}")


def plot():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams["svg.fonttype"] = "none"
    INK, MUTED, BLUE = "#1a1a1a", "#6b6b6b", "#0072b2"
    # Included at 0.49\textwidth of a 498.66pt (TeX) text block, the same slot
    # as panel B beside it, so point sizes below are the sizes on the page.
    W = 0.49 * 498.66 / 72.27

    d = pd.read_csv(TSV, sep="\t")
    assert len(d) == 475, f"{len(d)} rows, expected 475"
    x, y = d.gt_collisions.to_numpy(float), d.estimated_collisions.to_numpy(float)
    slope, intercept = np.polyfit(x, y, 1)
    r2 = np.corrcoef(x, y)[0, 1] ** 2

    fig, ax = plt.subplots(figsize=(W, 1.9), layout="constrained")
    top = max(x.max(), y.max()) * 1.04
    ax.plot([0, top], [0, top], color=MUTED, lw=0.8, ls=(0, (3, 2)), zorder=1)
    ax.scatter(x, y, s=9, color=BLUE, edgecolor="white", linewidth=0.4, zorder=2)
    ax.set_xlim(0, top)
    ax.set_ylim(0, top)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("true collisions", fontsize=8, color=INK)
    ax.set_ylabel("estimated collisions", fontsize=8, color=INK)
    ax.xaxis.set_major_formatter(matplotlib.ticker.FuncFormatter(
        lambda v, _: f"{v / 1000:g}k" if v else "0"))
    ax.yaxis.set_major_formatter(ax.xaxis.get_major_formatter())
    ax.tick_params(labelsize=7, colors=INK, length=2)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    # Below the diagonal, the one region the points never reach. The dashed
    # identity line and the number of simulations are in the caption.
    ax.text(0.97, 0.04, f"slope {slope:.3f}\n$R^2$ = {r2:.5f}",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=7,
            color=INK, linespacing=1.3)
    fig.savefig(SVG, format="svg")
    plt.close(fig)
    print(f"slope {slope:.4f} intercept {intercept:.2f} R2 {r2:.5f}; wrote {SVG}")


if __name__ == "__main__":
    {"simulate": simulate, "plot": plot}[sys.argv[1]]()
