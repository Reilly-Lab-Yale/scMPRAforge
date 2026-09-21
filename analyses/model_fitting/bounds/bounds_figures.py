"""Shared canvas and type settings for the Bounds transfection/library plots.

These panels are included in the manuscript at 0.30\\textwidth (149.6pt of a
498.66pt text block). LaTeX scales a figure by (include width / canvas width)
and every label scales with it, so the canvas is drawn at its include width
here and the point sizes are set to the house convention: 7pt annotations,
8pt axis labels, 9pt panel titles.
"""
import re

import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

# 0.30 * 498.66pt, in inches at matplotlib's 72pt/in SVG unit.
PANEL_W_IN = 149.60 / 72
PANEL_H_IN = 2.1

RC = {
    "svg.fonttype": "none",    # keep <text> editable for downstream relabelling
    "font.size": 7,
    "axes.labelsize": 8,
    "axes.titlesize": 9,
    "xtick.labelsize": 7,
    "ytick.labelsize": 7,
    "legend.fontsize": 7,
}


def save_model_plots(bounds, short, out_dir):
    """Write <short>_transfection.svg and <short>_library.svg into out_dir.

    The title names the dataset and the model component only. A full ortho
    name ("cohen_obsingle_nb_phantom") overruns the panel width at 9pt, more
    so after downstream relabelling (cohen -> Zhao et al.); the model variant
    stays in the file name.
    """
    written = []
    dataset = short.split("_")[0]
    with plt.rc_context(RC):
        for kind, draw in (("transfection", bounds.plot_transfection),
                           ("library", bounds.library_model.plot)):
            fig, ax = plt.subplots(figsize=(PANEL_W_IN, PANEL_H_IN),
                                   layout="constrained")
            fig.get_layout_engine().set(w_pad=0.01, h_pad=0.01)
            draw(ax=ax)
            ax.set_title(f"{dataset}\n{kind} model")
            # Counts run to six digits, which collides at the default tick
            # density once the axis is only ~1.5in wide. Fewer ticks, and a
            # shared power-of-ten offset once the labels get long.
            ax.xaxis.set_major_locator(MaxNLocator(nbins=4))
            ax.ticklabel_format(axis="x", style="sci", scilimits=(0, 3))
            # The legend labels are nearly as wide as the panel, so they go
            # under the axes rather than over the data. A constrained layout
            # reserves the room; tight_layout does not see an outside legend.
            handles, labels = ax.get_legend_handles_labels()
            ax.get_legend().remove()
            fig.legend(handles, labels, loc="outside lower center",
                       frameon=False, handlelength=1.2, handletextpad=0.4,
                       borderpad=0.1, labelspacing=0.2)
            path = out_dir / f"{short}_{kind}.svg"
            fig.savefig(str(path))
            plt.close(fig)
            saved_pt = float(re.search(r'<svg[^>]*width="([\d.]+)pt"',
                                       path.read_text()).group(1))
            want_pt = PANEL_W_IN * 72
            assert abs(saved_pt - want_pt) < 0.01 * want_pt, (
                f"{path.name}: saved canvas {saved_pt:.1f}pt, "
                f"include width {want_pt:.1f}pt -- type would not render at "
                f"the size it was set to")
            written.append(path)
    return written
