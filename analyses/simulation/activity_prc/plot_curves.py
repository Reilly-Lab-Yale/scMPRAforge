#!/usr/bin/env python
"""Median-replicate ROC and PRC curves, one pair per dataset.

Scriptified from the paper-figure cell of all_prc_summary.ipynb so the curves
can be regenerated without stepping through a notebook. Two changes:

  - Wald is dropped from Yin et al. Wald is the only test that reads the
    fitted model rather than the counts, and those orthos were fit under plain
    consider-missing rather than with the MOI correction that dataset's
    canonical fit uses (design fit_mode='cm_phantom'), so its curve there is
    not comparable to the other regimes. Pseudobulk is absent from Yin et al.
    on its own: it needs at least two replicates per group and that design has
    one.
  - The curve shown is chosen the same way as before -- the ground-truth draw
    whose mean auROC is closest to the median across draws -- but the choice
    is asserted rather than assumed.

Needs the simulation objects, hence a Dask client; run it on the cluster.

    python analyses/simulation/activity_prc/plot_curves.py
"""
import pathlib
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.collections import PathCollection
import pandas as pd

sys.path.insert(0, "/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa")
import scMPRAforge as scm  # noqa: E402

plt.rcParams["svg.fonttype"] = "none"

BASE = pathlib.Path(__file__).resolve().parent
OUT = BASE / "output" / "paper_figs"
SIM_ROOT = pathlib.Path("/nfs/roberts/project/pi_skr2/shared/tabula_data_new/simulated")

N_GT_DRAWS = 5
HYPOTHESIS_SET = "hs_all_ct"

# The panel is drawn at exactly the width Fig. 2 includes it at, so LaTeX
# applies no scaling and the type sizes set here are the sizes on the page:
# 7pt for tick labels, legend and annotation, 8pt for axis labels, 9pt for the
# panel title. Margins are in inches, hence the explicit subplots_adjust --
# a tight bbox would crop the canvas back off its target width.
PANEL_W, PANEL_H = 3.04, 2.70
MARGIN_L, MARGIN_R, MARGIN_T, MARGIN_B = 0.50, 0.08, 0.26, 0.42
TICK_PT, LEGEND_PT, ANNOT_PT, AXIS_PT, TITLE_PT = 7, 7, 7, 8, 9

# Legends carry the internal test identifiers; the manuscript uses these.
TEST_LABEL = {"mwu": "MWU", "ttest": "t-test", "ks": "KS",
              "pseudobulk": "pseudobulk", "wald_auto": "Wald"}

# One legend serves several panels once these are assembled (Fig. 2A pairs ROC
# with PRC, Fig. S2 puts one legend over five panels), so it may carry only
# what every panel shares. AUC belongs to a ROC panel and AP to a PRC panel,
# and each panel has its own baseline, so all of that is stripped here and the
# baselines are labelled directly on the line they describe instead.
BASELINE_LABEL = re.compile(r"^(chance|baseline=([0-9.]+))$")

# Each curve type leaves a different corner empty.
LEGEND_LOC = {"ROC": "lower right", "PRC": "lower left"}

# A fixed hue per test, for the same reason the legend carries no numbers: it
# is shared. Matplotlib colours by draw order, so the Yin et al. panels, which
# run three tests rather than five, would put t-test in the hue the legend has
# already given to pseudobulk.
#
# Okabe-Ito, as in activity_calibration/fpr_dumbbell.py, which gives MWU and
# the t-test these same two hexes. Of the eight Okabe-Ito slots only blue,
# vermillion and bluish green clear the categorical checks against a white page
# unaltered; the rest sit under 3:1 contrast or outside the lightness band. The
# green is therefore re-stepped (a change no eye resolves) and the last two
# slots are darker steps on the purple side, where every Okabe-Ito entry is too
# light to clear contrast against the page.
#
# All five clear the dataviz checks on every pair, not just neighbours, which
# is what five curves sharing one axis need. MWU holds the widest margin over
# those thresholds of the five: it is the test the paper adopts and the curve a
# reader follows.
#
#   node scripts/validate_palette.js \
#       "#0072b2,#d55e00,#1ca271,#5d06ca,#90026f" --pairs all --mode light
TEST_COLOR = {"mwu": "#0072b2", "ttest": "#d55e00", "ks": "#1ca271",
              "pseudobulk": "#5d06ca", "wald_auto": "#90026f"}

# Chance and the PRC prevalence line are references, not results. Left to the
# property cycle they take the next colour after the curves, which on a
# five-test panel is a sixth hue nobody chose and on the others is a hue a
# curve already holds.
BASELINE_GREY = "#444444"

DATASETS = {
    "shendure": ("Lalanne et al.", "shendure_5x5_activity"),
    "cohen":    ("Zhao et al.",    "cohen_5x5_activity"),
    "seelig":   ("Yin et al.",     "seelig_5x5_activity"),
}
# Tests to draw per dataset. None means "whatever ran".
#
# Lalanne et al. and Zhao et al. both run all five. Neither is special, and the
# manuscript features Lalanne et al. in Fig. 2A only as a representative
# example -- Zhao et al. would do as well. Yin et al. is the one that draws
# three, and the two missing tests are missing for unrelated reasons:
#
#   pseudobulk  never ran. It collapses each element's cells to one value per
#               replicate and compares those, so it needs at least two
#               replicates per group, and that design has one. Absent from the
#               summary TSV entirely.
#   wald_auto   ran, and is dropped here on purpose. It is the only test that
#               reads the fitted model rather than the counts, and the orthos
#               refit for this benchmark used plain consider-missing rather
#               than the MOI correction the canonical fit uses
#               (fit_mode='cm_phantom'), so its variance estimates are not
#               comparable to the other two regimes. Its median auROC there is
#               0.633 against 0.83-0.84 for the rest, which is an artefact of
#               that mismatch and not a result about the test.
#
# plot_test_comparison.py drops the same two cells via EXCLUDED; keep the two
# in step.
TESTS_FOR = {
    "shendure": None,
    "cohen": None,
    "seelig": ["mwu", "ttest", "ks"],
}


def median_gt_draw(slug):
    """The draw whose mean auROC is nearest the median across draws."""
    f = BASE / slug / "output" / f"{slug}_5x5_activity_summary.tsv"
    assert f.is_file(), f"missing summary: {f}"
    d = pd.read_csv(f, sep="\t")
    gt_means = d.groupby("gt_draw")["auroc"].mean()
    assert len(gt_means) == N_GT_DRAWS, (
        f"{slug}: expected {N_GT_DRAWS} ground-truth draws, found {len(gt_means)}")
    pick = int((gt_means - gt_means.median()).abs().idxmin())
    print(f"  {slug}: GT draw {pick} "
          f"(mean auROC {gt_means[pick]:.3f}; draws span "
          f"{gt_means.min():.3f}-{gt_means.max():.3f})")
    return pick


def relabel(entry):
    """'wald_auto (AUC=0.929)' -> 'Wald'; other entries untouched."""
    head, _, _rest = entry.partition(" ")
    return TEST_LABEL.get(head, entry)


def recolor(ax):
    """Give each test its own hue regardless of how many ran on this panel."""
    curves, refs = [], []
    for ln in ax.get_lines():
        (curves if ln.get_label().partition(" ")[0] in TEST_COLOR else refs).append(ln)
    assert curves, f"no test curves found among {[l.get_label() for l in ax.get_lines()]}"
    # One scatter per curve, added in curve order. The PRC baseline is a
    # LineCollection and must not be counted among them.
    marks = [c for c in ax.collections if isinstance(c, PathCollection)]
    assert len(marks) in (0, len(curves)), \
        f"{len(marks)} threshold markers against {len(curves)} curves"
    for i, line in enumerate(curves):
        colour = TEST_COLOR[line.get_label().partition(" ")[0]]
        line.set_color(colour)
        if marks:
            marks[i].set_color(colour)
    for ref in refs:
        ref.set_color(BASELINE_GREY)
    for coll in ax.collections:
        if not isinstance(coll, PathCollection):
            coll.set_color(BASELINE_GREY)
    drawn = {c.get_color() for c in curves}
    assert len(drawn) == len(curves), f"two curves share a hue: {drawn}"
    assert BASELINE_GREY not in drawn, "a curve is wearing the baseline grey"


def label_baseline(ax, kind, value):
    """Name the baseline on the line itself, where its panel can be seen."""
    grey = BASELINE_GREY
    if kind == "ROC":
        # Match the diagonal's on-screen slope; the axes is not square, so 45
        # degrees would visibly miss it.
        (x0, y0), (x1, y1) = ax.transData.transform([[0, 0], [1, 1]])
        ax.text(0.63, 0.60, "chance", fontsize=ANNOT_PT, color=grey,
                rotation=np.degrees(np.arctan2(y1 - y0, x1 - x0)),
                rotation_mode="anchor", ha="center", va="top")
    else:
        # Mid-line: the curves come down onto the right end and the legend
        # takes the left.
        ax.annotate(f"baseline {value:.3f}", xy=(0.62, value),
                    xytext=(0, 3), textcoords="offset points",
                    fontsize=ANNOT_PT, color=grey, ha="center", va="bottom")


def style_panel(fig, label, kind, slug=""):
    """Size the canvas, restyle the type, recolour, and fix up the legend."""
    fig.set_size_inches(PANEL_W, PANEL_H)
    fig.set_layout_engine("none")
    fig.subplots_adjust(left=MARGIN_L / PANEL_W, right=1 - MARGIN_R / PANEL_W,
                        top=1 - MARGIN_T / PANEL_H, bottom=MARGIN_B / PANEL_H)
    for ax in fig.axes:
        # One title, naming the regime; the generic "ROC Curve" the plotting
        # call sets is redundant beside it.
        ax.set_title(f"{label}, {kind}", fontsize=TITLE_PT)
        ax.xaxis.label.set_fontsize(AXIS_PT)
        ax.yaxis.label.set_fontsize(AXIS_PT)
        ax.tick_params(labelsize=TICK_PT)
        recolor(ax)
        if ax.get_legend() is None:
            continue
        handles, labels = ax.get_legend_handles_labels()
        base = [BASELINE_LABEL.match(t) for t in labels]
        assert sum(m is not None for m in base) == 1, (
            f"{slug} {kind}: expected exactly one baseline entry in {labels}")
        hit = next(m for m in base if m)
        label_baseline(ax, kind, float(hit.group(2) or 0.0))

        keep = [(h, relabel(t)) for h, t, m in zip(handles, labels, base)
                if m is None]
        assert keep, f"{slug} {kind}: every legend entry was dropped"
        ax.legend(*zip(*keep), loc=LEGEND_LOC[kind], fontsize=LEGEND_PT,
                  framealpha=0.95, borderpad=0.5, labelspacing=0.35,
                  handlelength=1.4)


def main():
    from dask.distributed import Client, LocalCluster
    cluster = LocalCluster(n_workers=1, threads_per_worker=1, memory_limit="16GB")
    client = Client(cluster)
    OUT.mkdir(parents=True, exist_ok=True)

    for slug, (label, sim_stem) in DATASETS.items():
        gt = median_gt_draw(slug)
        sim = scm.de_novo_simulation(
            location=SIM_ROOT, name=f"{sim_stem}_gt{gt}", client=client)
        tests = TESTS_FOR[slug]
        for kind in ("ROC", "PRC"):
            sim.median_performance_curve(
                HYPOTHESIS_SET, kind, test_types=tests, include_alpha=True)
            fig = plt.gcf()
            style_panel(fig, label, kind, slug)
            out = OUT / f"{slug}_median_{kind.lower()}.svg"
            fig.savefig(out, format="svg")
            plt.close(fig)
            print(f"    wrote {out.name}")

    client.close()
    print("done")


if __name__ == "__main__":
    main()
