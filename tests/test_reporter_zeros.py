"""Reporter-informed zero counts, per (cell, CRE) case.

A coarse reporter detects the CRE, not the barcode, so the two expansions
answer different questions and a (cell, CRE) pair falls into one of four
cases. The rule each expansion applies to those cases is the whole content
of `_reporter_zero_counts`, and getting one case wrong is silent: the fit
converges either way, only the fitted means move.

    python tests/test_reporter_zeros.py     (or: pytest tests/)
"""
import pandas as pd

import scMPRAforge.core as core

# CRE A carries three barcodes, CRE B two; one cell type, one replicate.
MPRA_MAP = pd.DataFrame({"rep_id": ["1"] * 5,
                         "mpra_bc": ["a1", "a2", "a3", "b1", "b2"],
                         "cre_id": ["A", "A", "A", "B", "B"]})
CELL_MAP = pd.DataFrame({"rep_id": ["1"] * 4,
                         "cell_bc": ["c1", "c2", "c3", "c4"],
                         "cell_type": ["T"] * 4})
# The four cases, one per cell:
#   c1  +reporter -MPRA   silent, the only case that earns a zero
#   c2  +reporter +MPRA   an expression measurement, not a silent one
#   c3  -reporter +MPRA   orphan; its own counts prove it was transfected
#   c4  +reporter -MPRA   silent, on CRE B
NONZERO = pd.DataFrame({"rep_id": ["1", "1"], "cell_bc": ["c2", "c3"],
                        "cre_id": ["A", "A"], "cell_type": ["T", "T"],
                        "mpra_umis": [4, 7]})
REPORTER = pd.DataFrame({"rep_id": ["1", "1", "1"],
                         "cell_bc": ["c1", "c2", "c4"],
                         "cre_id": ["A", "A", "B"]})


def n_total(expansion):
    t = core._reporter_zero_counts(NONZERO, REPORTER, MPRA_MAP, CELL_MAP,
                                   "cre_id", levels=["A", "B"],
                                   reporter_expansion=expansion)
    return dict(zip(t["cre_id"], t["n_total"]))


def test_single_counts_only_silent_detections():
    """One zero per detected-but-silent pair, which is Zhao et al.'s U.

    Here that is c1 on A and c4 on B. Counting every pair with evidence of
    transfection instead would give A three zeros and B one, inflating the
    zeros by one for each pair that expressed.
    """
    assert n_total("single") == {"A": 1, "B": 1}


def test_coarse_counts_candidates_including_orphans():
    """Candidate observations, from which nonzeros are subtracted downstream.

    Every pair with evidence of transfection contributes its CRE's barcodes,
    so A spans c1, c2 and c3 at three barcodes each. The orphan c3 belongs
    here: the reporter missed it, but its counts prove delivery.
    """
    assert n_total("coarse") == {"A": 9, "B": 2}




def test_single_keeps_groups_whose_pairs_all_expressed():
    """A group with no silent pair still gets a row, at zero.

    The caller reads its group list from this table and merges the nonzero
    counts onto it, so a dropped group would silently take its own
    observations out of the design. CRE C is detected once, in a cell that
    expressed it, so it earns no zero but must still appear.
    """
    reporter = pd.concat([REPORTER, pd.DataFrame(
        {"rep_id": ["1"], "cell_bc": ["c3"], "cre_id": ["C"]})])
    mpra_map = pd.concat([MPRA_MAP, pd.DataFrame(
        {"rep_id": ["1"], "mpra_bc": ["z1"], "cre_id": ["C"]})])
    nonzero = pd.concat([NONZERO, pd.DataFrame(
        {"rep_id": ["1"], "cell_bc": ["c3"], "cre_id": ["C"],
         "cell_type": ["T"], "mpra_umis": [2]})])
    t = core._reporter_zero_counts(nonzero, reporter, mpra_map, CELL_MAP,
                                   "cre_id", levels=["A", "B", "C"],
                                   reporter_expansion="single")
    got = dict(zip(t["cre_id"], t["n_total"]))
    assert got == {"A": 1, "B": 1, "C": 0}, got


if __name__ == "__main__":
    for name, fn in sorted(globals().items()):
        if name.startswith("test_"):
            fn()
            print(f"ok  {name}")
