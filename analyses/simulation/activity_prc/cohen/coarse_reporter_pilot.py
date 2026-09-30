"""
Coarse-reporter pilot: does collapsing simulated zeros to what a CRE-level
transfection reporter can observe change MWU performance in the Zhao et al.
regime?

The simulator emits one observation per (cell, delivered MPRA barcode), and
with a reporter every zero is kept. Zhao et al.'s U6 reporter only says which
CRE reached a cell, not which barcode, so a real assay cannot see those
per-barcode zeros. In this regime a delivered (cell, CRE) pair carries ~4.5
distinct barcodes on average, so the difference is not small in principle.

Re-runs MWU on the cached 5x5 activity simulations under four views of the
same counts, grouped by (rep_id, cell_bc, cell_type, cre_id):

  mwu       as simulated: every zero kept (the published arm)
  strict    nonzero observations, plus one zero if the pair has none
  cohen     nonzero observations, plus one zero for every delivered pair
            (Zhao et al.'s convention, and the canonical fit's "single")
  deflated  zeros dropped (no reporter)

Nothing in the package changes: collapse happens in a bundle_fn passed to
the package's own _counts_test_worker, and results are written here, not
into the shared simulation directories.

The mwu arm must reproduce the stored summary exactly; that checks the
scoring path before the other arms are read.

    python coarse_reporter_pilot.py
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score

import scMPRAforge.core as core

SIM_ROOT = Path("/nfs/roberts/project/pi_skr2/shared/tabula_data_new/simulated")
PREFIX = "cohen_5x5_activity"
HS = "hs_all_ct"
N_GT, N_REP = 5, 5
ALPHA = 0.05
HERE = Path(__file__).resolve().parent
OUT = HERE / "output" / "coarse_reporter_pilot"
STORED = HERE / "output" / "cohen_5x5_activity_summary.tsv"
PAIR = ["rep_id", "cell_bc", "cell_type", "cre_id"]
COLS = PAIR + ["umis_mpra_bc"]


def collapse(pdf, zero_rule):
    """Keep nonzero rows; add one zero per pair ('strict': only pairs with
    no nonzero row, 'cohen': every pair)."""
    n_pairs = len(pdf[PAIR].drop_duplicates())
    nz = pdf[pdf["umis_mpra_bc"] > 0]
    pairs = pdf[PAIR].drop_duplicates()
    if zero_rule == "strict":
        seen = nz[PAIR].drop_duplicates().assign(_nz=True)
        pairs = pairs.merge(seen, on=PAIR, how="left")
        pairs = pairs[pairs["_nz"].isna()].drop(columns="_nz")
    zeros = pairs.assign(umis_mpra_bc=0)
    out = pd.concat([nz[COLS], zeros[COLS]], ignore_index=True)
    # every delivered pair survives, and no pair gains more than one zero
    assert len(out[PAIR].drop_duplicates()) == n_pairs, \
        f"collapse lost pairs: {n_pairs} -> {len(out[PAIR].drop_duplicates())}"
    assert (out["umis_mpra_bc"] == 0).sum() <= n_pairs
    return out


def bundle_strict(pdf):
    return core._build_mwu_counts_dict(collapse(pdf, "strict"))


def bundle_cohen(pdf):
    return core._build_mwu_counts_dict(collapse(pdf, "cohen"))


ARMS = {  # arm -> (has_reporter, bundle_fn)
    "mwu": (True, None),
    "strict": (True, bundle_strict),
    "cohen": (True, bundle_cohen),
    "deflated": (False, None),
}


def score(res, gt):
    """AUROC/AUPRC exactly as de_novo_simulation._classifier_summary, plus
    the false-positive rate on true nulls at BH < ALPHA."""
    gt = gt[["cre_id", "cell_type", "true_mean"]]
    m = res.merge(gt, left_on=["comparison_CRE", "comparison_cell_type"],
                  right_on=["cre_id", "cell_type"]).drop(columns=["cre_id", "cell_type"])
    m = m.rename(columns={"true_mean": "comparison_truth"})
    m = m.merge(gt, left_on=["reference_CRE", "reference_cell_type"],
                right_on=["cre_id", "cell_type"]).drop(columns=["cre_id", "cell_type"])
    m = m.rename(columns={"true_mean": "reference_truth"})
    assert len(m) == len(res), f"ground-truth merge changed rows: {len(res)} -> {len(m)}"
    m["gt_null"] = (m["comparison_truth"] / m["reference_truth"] - 1).abs() < core.FLOATING_POINT_DIFF
    m["meta"] = m["meta"].astype("string").fillna("0")
    m = m.dropna()
    y = (~m["gt_null"]).astype(int)
    s = 1.0 - m["p_value"]
    null = m[m["gt_null"]]
    return {"auroc": roc_auc_score(y, s), "auprc": average_precision_score(y, s),
            "fpr": float((null["bh_p"] < ALPHA).mean()), "n_null": len(null),
            "n_tests": len(m)}


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    rows = []
    for g in range(N_GT):
        simd = SIM_ROOT / f"{PREFIX}_gt{g}"
        hyp = core.HypothesisSet.from_tsv(simd / "tests" / HS / "hypotheses.tsv")
        gt = pd.read_csv(simd / "ground_truth.tsv.gz", sep="\t", index_col=0)
        gt = core.cast_string_keys(gt.rename(columns={"mu": "true_mean"}), ["cell_type", "cre_id"])
        for r in range(N_REP):
            for arm, (has_rep, bfn) in ARMS.items():
                f = OUT / f"gt{g}_{r}_{arm}.tsv"
                if not f.is_file():
                    core._counts_test_worker(
                        None, simd / "simulated_scmpra" / f"{r}.scmpra", f, hyp,
                        has_rep, core._mwu_row_fn, f"mwu_{arm}",
                        columns=COLS, bundle_fn=bfn)
                res = core.ResultSet.from_tsv(f).df
                rows.append({"gt_draw": g, "replicate": r, "arm": arm, **score(res, gt)})
            print(f"gt{g} rep{r} done", flush=True)

    out = pd.DataFrame(rows)
    out.to_csv(OUT / "summary.tsv", sep="\t", index=False)

    # Scoring check: the as-simulated arm must match the published summary.
    stored = pd.read_csv(STORED, sep="\t").query("test == 'mwu'")
    chk = out.query("arm == 'mwu'").merge(stored, on=["gt_draw", "replicate"],
                                          suffixes=("", "_stored"))
    assert len(chk) == N_GT * N_REP, f"stored summary match: {len(chk)} rows"
    d = (chk["auprc"] - chk["auprc_stored"]).abs().max()
    assert d < 1e-9, f"scoring does not reproduce stored mwu auPRC (max diff {d})"
    print(f"scoring check: mwu arm reproduces stored auPRC (max diff {d:.1e})")

    print(out.groupby("arm")[["auprc", "auroc", "fpr"]]
          .agg(["median", "min", "max"]).round(3).to_string())


if __name__ == "__main__":
    sys.exit(main())
