#!/usr/bin/env python
"""Estimated collision rate per dataset -- Fig S1B.

A collision is a cell receiving a barcode it already carries, which is what
the RNA-only model assumes away. The three numbers are the estimator's output
on each empirical dataset, not a fit, so they are held here as literals rather
than recomputed: rerunning the estimator needs the full UMI tables.

Provenance, matching the notebook this replaces:

    Zhao et al.     0.11689209661257542   INFO from Bounds.from_ortho,
    Lalanne et al.  0.05304031883524153   excess / total_tfection * 100
    Yin et al.      7.630123998189407     seelig_collision_rate.py, raw TSV

Yin et al. is two orders of magnitude higher because its barcodes are
synthesised rather than appended as random N-mers, so its library carries
about 5 barcodes per element against Lalanne et al.'s 136.

    python analyses/supplementary/plot_collision_rate.py
"""
import pathlib

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Editable <text> so the manuscript sync pipeline's transforms apply.
plt.rcParams["svg.fonttype"] = "none"

BASE = pathlib.Path(__file__).resolve().parent
OUT = BASE / "estimated_percent_conflict_collision_rate.svg"

# Included at 0.49\textwidth of a 498.66pt (TeX) text block. matplotlib works
# in PostScript points, so the check below is against that, not TeX points.
DISPLAY_W_IN = 0.49 * 498.66 / 72.27
BLUE, INK, MUTED = "#0072b2", "#1a1a1a", "#6b6b6b"   # Okabe-Ito, house palette

PERCENT_COLLISION = {
    "Zhao et al.": 0.11689209661257542,
    "Lalanne et al.": 0.05304031883524153,
    "Yin et al.": 7.630123998189407,
}


def main():
    rows = sorted(PERCENT_COLLISION.items(), key=lambda kv: -kv[1])
    names = [k for k, _ in rows]
    vals = [v for _, v in rows]
    assert all(v >= 0 for v in vals), f"negative collision rate: {rows}"
    assert max(vals) < 100, f"collision rate above 100%: {rows}"

    # Constrained layout rather than a tight bbox: a tight bbox trims or
    # expands the saved canvas away from figsize, and the whole point here is
    # that the saved width equals the width the figure is included at.
    fig, ax = plt.subplots(figsize=(DISPLAY_W_IN, 1.9), layout="constrained")
    ax.bar(names, vals, color=BLUE, width=0.62)

    # One series, so no legend -- the axis label names the quantity.
    ax.set_ylabel("collision rate\n(% of transfection events)", fontsize=8,
                  color=INK)
    ax.set_ylim(0, max(vals) * 1.22)
    ax.bar_label(ax.containers[0], labels=[f"{v:.3f}%" for v in vals],
                 fontsize=7, color=MUTED, padding=2)

    ax.grid(True, axis="y", color="#e6e6e6", lw=0.8, zorder=0)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color("#cccccc")
    ax.tick_params(colors=MUTED, labelsize=7, length=0)

    fig.savefig(OUT, format="svg")
    plt.close(fig)

    # The saved width is what LaTeX scales against; drift here silently
    # rescales every label on the page.
    import re
    saved = float(re.search(r'width="([\d.]+)pt"', OUT.read_text()).group(1))
    want = DISPLAY_W_IN * 72
    assert abs(saved - want) / want < 0.01, (
        f"saved {saved:.1f}pt but the slot is {want:.1f}pt; "
        "the type will not render at the size set here")
    print(f"wrote {OUT} ({saved:.1f}pt = {saved / 72:.2f} in)")


if __name__ == "__main__":
    main()
