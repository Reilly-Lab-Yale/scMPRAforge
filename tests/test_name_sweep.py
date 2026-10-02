"""Legacy file boundaries and canonical derived-count consumers."""
import importlib.util
from pathlib import Path

import pandas as pd
import pytest

from scMPRAforge import core

ROOT = Path(__file__).resolve().parents[1]


def load_script(relative):
    path = ROOT / relative
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize("name", ["mpra_umis", "umis_mpra_bc"])
def test_standalone_readers_accept_both_names(tmp_path, name):
    pdf = pd.DataFrame({"cre_id": ["a", "b"], "cell_type": ["ct", "ct"], name: [0, 5]})
    tsv = tmp_path / "counts.tsv"
    pdf.to_csv(tsv, sep="\t", index=False)
    overdispersion = load_script("analyses/model_selection/overdispersion.py")
    got = overdispersion.read_counts(tsv)
    assert got["mpra_umis"].tolist() == [0, 5]
    assert len(got) == len(pdf)
    backfill = load_script("analyses/model_fitting/backfill_ortho_meta.py")
    assert backfill.source_has_zeros(tsv) is True
    scmpra = tmp_path / "counts.scmpra"
    scmpra.mkdir()
    pdf.to_parquet(scmpra / "data.parquet", index=False)
    assert backfill.source_has_zeros(scmpra) is True


def test_standalone_readers_reject_competing_names(tmp_path):
    path = tmp_path / "counts.tsv"
    pd.DataFrame({"cre_id": ["a"], "cell_type": ["ct"],
                  "mpra_umis": [2], "umis_mpra_bc": [3]}).to_csv(path, sep="\t", index=False)
    overdispersion = load_script("analyses/model_selection/overdispersion.py")
    backfill = load_script("analyses/model_fitting/backfill_ortho_meta.py")
    with pytest.raises(AssertionError, match="one MPRA count column"):
        overdispersion.read_counts(path)
    with pytest.raises(AssertionError, match="one MPRA count column"):
        backfill.source_has_zeros(path)


@pytest.mark.parametrize("name", ["mpra_umis_normalized", "normalized_umis_mpra_bc"])
def test_bootstrap_uses_normalized_counts_and_preserves_them(tmp_path, name):
    data = core.scMPRA_data()
    data.data = pd.DataFrame({"rep_id": ["r", "r"], "cell_bc": ["c1", "c2"],
                              "cell_type": ["ct", "ct"], "cre_id": ["reference", "active"],
                              "transfection_bc": ["t1", "t2"], "mpra_umis": [10, 20],
                              name: [0.5, 1.25]})
    data.table_type = "mpra_umiwise"
    hs = core.HypothesisSet.from_dataframe(pd.DataFrame({
        "comparison_CRE": ["active"], "comparison_cell_type": ["ct"],
        "reference_CRE": ["reference"], "reference_cell_type": ["ct"],
    }))
    bundle = core._bootstrap_build_bundle(hs, data)
    assert bundle["metric_col"] == "mpra_umis_normalized"
    assert sorted(bundle["integrations"]["value"]) == [0.5, 1.25]
    dest = tmp_path / "normalized.scmpra"
    data.to_parquet(dest)
    restored = core.scMPRA_data.from_parquet(dest)
    assert restored.data.compute()["mpra_umis_normalized"].tolist() == [0.5, 1.25]
    assert "normalized_umis_mpra_bc" not in restored.data.columns


def test_coarse_pilot_passes_canonical_counts_to_worker_bundle():
    pilot = load_script("analyses/simulation/activity_prc/cohen/coarse_reporter_pilot.py")
    pdf = pd.DataFrame({"rep_id": ["r"] * 3, "cell_bc": ["c1", "c1", "c2"],
                        "cell_type": ["ct"] * 3, "cre_id": ["a"] * 3,
                        "mpra_umis": [2, 0, 0]})
    bundle = pilot.bundle_strict(pdf)
    assert sorted(bundle["counts"][("ct", "a")]) == [0, 2]
