"""Re-derive the EMPIRICAL anchors of the design-space sweep from the fits.

The anchors place each published design as a point in the swept space, and
every one of them is a number read off that dataset's canonical Bounds preset
and canonical ortho. They live in synthetic_factorial.py as literals, so a
refit silently leaves them describing the previous fit. This script derives
them again and diffs against what is written there.

Two of the seven come from the fitted activities and so move whenever the
canonical ortho is refit; the rest come from the library, transfection and
cell counts, which a refit of the same data leaves alone:

    minP               reference_activity                     (fit)
    activity_max_mult  p95(mu in the reference cell type)/minP (fit)
    n_cells            cells_per_cell_type['reference']
    n_cres             CREs in the reference cell type
    bcs_per_cre        library_model.mu_nb
    lib_alpha_nb       library_model.alpha_nb
    moi                transfection_model.mu_nb

Run it after any canonical refit. It exits 1 if anything drifted, and prints
a dict ready to paste over EMPIRICAL.

    python derive_anchors.py            (on Bouchet; needs the shared data)
"""
import pickle
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, "/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa")
import scMPRAforge.core as scm  # noqa: E402

from synthetic_factorial import EMPIRICAL  # noqa: E402

PRESETS = Path("/nfs/roberts/project/pi_skr2/mcn26/tabula-rasa/scMPRAforge/presets")
DATA = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data_new")

# anchor name -> (dataset, canonical fit). The fit named here is the one the
# manuscript reports for that dataset; an anchor derived from any other fit
# would place the design somewhere the paper never claims it is.
CANONICAL = {
    "shendure-Pluripotent": ("shendure", "shendure_obs_nb"),
    "cohen-Rod": ("cohen", "cohen_obsingle_nb_phantom"),
    "seelig-HepG2": ("seelig", "seelig_cm_moib_nb_phantom"),
}

# Drift below this is rounding in the literals, not a changed fit.
TOL = 5e-3


def derive(dataset, fit):
    b = scm.Bounds.from_tgz(str(PRESETS / f"{fit}.tgz"))
    params = pickle.load(
        open(DATA / dataset / fit / "by_cell_type_parameters.pkl", "rb"))
    mu = np.asarray(params.nb["reference"]["mu"])
    assert mu.size and np.isfinite(mu).all(), f"{fit}: bad mu in reference"
    ref = float(b.reference_activity)
    assert ref > 0, f"{fit}: reference_activity is {ref}"
    # The preset and the ortho are separate artifacts and a refit updates them
    # in two steps, so a stale preset read against a fresh ortho would yield a
    # baseline from one fit and a spread from the other. reference_activity is
    # the reference CRE's mu averaged over cell types, which ties the two
    # together and fails loudly when only one has been regenerated.
    across = np.array([float(params.nb[k]["mu"]["reference"]) for k in params.keys])
    assert abs(across.mean() - ref) / ref < 1e-3, (
        f"{fit}: preset reference_activity {ref:.6g} does not match the ortho's "
        f"mean reference mu {across.mean():.6g}; one of them is stale -- "
        f"regenerate the Bounds preset from this ortho before deriving anchors")
    return dict(
        n_cells=int(b.cells_per_cell_type["reference"]),
        n_cres=int(mu.size),
        bcs_per_cre=float(b.library_model.mu_nb),
        moi=float(b.transfection_model.mu_nb),
        lib_alpha_nb=float(b.library_model.alpha_nb),
        minP=ref,
        activity_max_mult=float(np.percentile(mu, 95) / ref),
    )


def main():
    drift = []
    out = {}
    for anchor, (dataset, fit) in CANONICAL.items():
        got = derive(dataset, fit)
        out[anchor] = got
        have = EMPIRICAL.get(anchor, {})
        print(f"\n== {anchor}  ({fit})")
        for k, v in got.items():
            old = have.get(k)
            if old is None:
                print(f"   {k:18s} {v:>12.4g}   (not in EMPIRICAL)")
                continue
            rel = abs(v - old) / max(abs(old), 1e-12)
            flag = "" if rel <= TOL else "   <-- DRIFT"
            print(f"   {k:18s} {v:>12.4g}   was {old:<12.4g}{flag}")
            if rel > TOL:
                drift.append(f"{anchor}.{k}: {old} -> {v:.6g}")

    if drift:
        print("\nEMPIRICAL is out of date:")
        for d in drift:
            print("  " + d)
        print("\nPaste over EMPIRICAL in synthetic_factorial.py:\n")
        for anchor, got in out.items():
            body = ", ".join(f"{k}={v:.6g}" if isinstance(v, float) else f"{k}={v}"
                             for k, v in got.items())
            print(f'    "{anchor}": dict({body}),')
        return 1
    print("\nEMPIRICAL matches the canonical fits.")
    return 0


if __name__ == "__main__":
    sys.exit(main())
