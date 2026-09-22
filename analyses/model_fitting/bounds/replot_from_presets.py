"""Redraw the manuscript transfection/library panels from the saved presets.

The extract_*_bounds.py scripts fit the orthos and save both the presets and
the plots; this redraws the plots alone, so a figure change does not need the
fits re-run. The preset round-trips the fitted models, so the panels are the
same plot.

    python analyses/model_fitting/bounds/replot_from_presets.py
"""
import matplotlib
matplotlib.use("Agg")

from pathlib import Path

import bounds_figures
import scMPRAforge.core as scm

HERE = Path(__file__).resolve().parent
PRESET_DIR = HERE.parents[2] / "scMPRAforge" / "presets"

# The canonical ortho per dataset (paper_plan.md), and the output subdirectory
# the extract script writes it to.
PANELS = [
    ("shendure", "shendure_obs_nb"),
    ("cohen", "cohen_obsingle_nb_phantom"),
    ("seelig", "seelig_cm_moib_nb_phantom"),
]


def main():
    for subdir, short in PANELS:
        out_dir = HERE / "output" / subdir
        out_dir.mkdir(parents=True, exist_ok=True)
        bounds = scm.Bounds.from_tgz(str(PRESET_DIR / f"{short}.tgz"))
        for path in bounds_figures.save_model_plots(bounds, short, out_dir):
            print(f"Saved: {path}", flush=True)


if __name__ == "__main__":
    main()
