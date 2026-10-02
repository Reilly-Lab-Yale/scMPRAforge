"""Canonical schemas, legacy imports, and count-preserving conversion."""
import json
import pickle

import dask.dataframe as dd
import pandas as pd
import pytest

from scMPRAforge import core


def readwise():
    return pd.DataFrame({
        "rep_id": ["01"] * 3,
        "cell_bc": ["001", "001", "002"],
        "cell_type": ["NA"] * 3,
        "cre_id": ["cre1"] * 3,
        "mpra_bc": ["bc1"] * 3,
        "mpra_umi": ["u1", "u2", "u1"],
        "mpra_reads": [2, 3, 4],
        "transfection_bc": ["t1", "t1", "t2"],
        "transfection_umis": [2, 2, None],
        "transfection_reads": [5, 5, None],
        "dna_reads": [20, 20, 20],
    })


def obj(pdf, kind="mpra_readwise"):
    data = core.scMPRA_data()
    data.table_type = kind
    data.data = dd.from_pandas(pdf, npartitions=2)
    return data


def legacy(pdf):
    return pdf.rename(columns={"mpra_umi": "umi", "mpra_reads": "reads",
                               "mpra_umis": "umis_mpra_bc",
                               "transfection_umis": "umis_transfection_bc",
                               "transfection_reads": "reads_transfection_bc",
                               "dna_reads": "reads_DNA"})


@pytest.mark.parametrize("old_names", [False, True])
def test_tsv_conversion_and_parquet_roundtrip(tmp_path, old_names):
    pdf = readwise()
    path = tmp_path / "reads.tsv"
    (legacy(pdf) if old_names else pdf).to_csv(path, sep="\t", index=False)
    data = core.scMPRA_data.from_tsv(path)
    assert set(data.data.columns) == set(pdf.columns)
    data.read_wise_to_umi_wise(keep_reads=True)
    out = data.data.compute().sort_values("cell_bc").reset_index(drop=True)
    assert out["rep_id"].tolist() == ["01", "01"]
    assert out["cell_bc"].tolist() == ["001", "002"]
    assert out["cell_type"].tolist() == ["NA", "NA"]
    assert out["mpra_umis"].sparse.to_dense().tolist() == [2, 1]
    assert out["mpra_reads"].tolist() == [5, 4]
    assert out.loc[0, "transfection_umis"] == 2
    assert out.loc[0, "transfection_reads"] == 5
    assert pd.isna(out.loc[1, "transfection_umis"])
    assert pd.isna(out.loc[1, "transfection_reads"])
    assert out["dna_reads"].tolist() == [20, 20]
    data.set_coarse_reporter(pd.DataFrame({"rep_id": ["01", "01"],
                                         "cell_bc": ["001", "003"], "cre_id": ["cre1", "cre1"]}))
    dest = tmp_path / "data.scmpra"
    data.to_parquet(dest)
    members = json.loads((dest / "members.json").read_text())
    assert members["schema_version"] == 2
    assert (dest / "members.pkl").exists()
    disk = pd.read_parquet(dest / "data.parquet")
    assert not set(disk.columns).intersection(core.MPRA_LEGACY_COLUMNS)
    loaded = core.scMPRA_data.from_parquet(dest).validate()
    restored = loaded.data.compute().sort_values("cell_bc").reset_index(drop=True)
    pd.testing.assert_frame_equal(out, restored[out.columns])
    pd.testing.assert_frame_equal(data._coarse_reporter, loaded._coarse_reporter)


def test_readwise_conversion_without_reads_and_separate_integrations():
    pdf = readwise()
    pdf.loc[1, "mpra_umi"] = "u1"
    pdf.loc[1, "transfection_bc"] = "different_integration"
    data = obj(pdf)
    data.read_wise_to_umi_wise()
    out = data.data.compute()
    assert "mpra_reads" not in out
    assert "mpra_umi" not in out
    assert len(out) == 3
    assert out["mpra_umis"].sum() == 3


@pytest.mark.parametrize("column, value", [
    ("mpra_reads", None), ("mpra_reads", 0), ("mpra_reads", -1),
    ("mpra_reads", 1.5), ("mpra_reads", float("inf")), ("mpra_reads", "bad"),
    ("transfection_umis", -1), ("transfection_reads", 1),
    ("cell_bc", None), ("mpra_umi", ""),
])
def test_invalid_readwise_values_fail(column, value):
    pdf = readwise().astype({column: object})
    pdf.loc[0, column] = value
    with pytest.raises((AssertionError, ValueError, TypeError)):
        obj(pdf).validate()


@pytest.mark.parametrize("column, value", [
    ("transfection_umis", 3), ("transfection_reads", 6),
    ("dna_reads", 21), ("cell_type", "other"), ("cre_id", "other"),
])
def test_conflicting_annotations_fail_across_partitions(column, value):
    pdf = readwise()
    pdf.loc[1, column] = value
    with pytest.raises(AssertionError, match="conflicting"):
        obj(pdf).read_wise_to_umi_wise(True)


def test_duplicate_umi_is_rejected():
    pdf = readwise()
    pdf.loc[1, "mpra_umi"] = "u1"
    with pytest.raises(AssertionError, match="duplicate observation keys"):
        obj(pdf).read_wise_to_umi_wise()


@pytest.mark.parametrize("columns", [
    {"reads": [2, 3, 4]}, {"reads_mpra_bc": [2, 3, 4]},
    {"transfection_umi": ["x", "y", "z"]},
])
def test_ambiguous_or_unpaired_legacy_columns_rejected(tmp_path, columns):
    pdf = readwise().assign(**columns)
    path = tmp_path / "bad.tsv"
    pdf.to_csv(path, sep="\t", index=False)
    with pytest.raises(ValueError):
        core.scMPRA_data.from_tsv(path)


def test_legacy_parquet_projection_and_pickle(tmp_path):
    pdf = pd.DataFrame({"rep_id": ["r"], "cell_type": ["t"],
                        "cre_id": ["c"], "umis_mpra_bc": [2], "reads_mpra_bc": [4]})
    dest = tmp_path / "legacy.scmpra"
    dest.mkdir()
    pdf.to_parquet(dest / "data.parquet", index=False)
    (dest / "members.json").write_text(json.dumps({"table_type": "mpra_umiwise"}))
    data = core.scMPRA_data.from_parquet(dest).validate()
    assert data.data.compute()["mpra_umis"].iloc[0] == 2
    projected = core._read_mpra_parquet(dest / "data.parquet", ["cell_type", "cre_id", "mpra_umis"])
    assert projected.columns.tolist() == ["cell_type", "cre_id", "mpra_umis"]
    assert projected["mpra_umis"].tolist() == [2]
    old = obj(pdf, "mpra_umiwise")
    restored = pickle.loads(pickle.dumps(old))
    assert "mpra_umis" in restored.data.columns
    assert "umis_mpra_bc" not in restored.data.columns


@pytest.mark.parametrize("bad", [None, -1, 1.5, "bad"])
def test_required_umi_counts_are_not_coerced_to_zero(bad):
    pdf = pd.DataFrame({"rep_id": ["r"], "cell_type": ["t"], "cre_id": ["c"], "mpra_umis": [bad]})
    with pytest.raises((ValueError, AssertionError)):
        obj(pdf, "mpra_umiwise").validate()


def test_anonymous_simulation_rows_are_not_deduplicated():
    pdf = pd.DataFrame({"rep_id": ["r"] * 2, "cell_type": ["t"] * 2,
                        "cre_id": ["c"] * 2, "mpra_umis": [0, 0]})
    data = obj(pdf, "mpra_umiwise").validate()
    assert len(data.data.compute()) == 2


def test_dna_counts_require_barcode():
    pdf = pd.DataFrame({"rep_id": ["r"], "cell_type": ["t"],
                        "cre_id": ["c"], "mpra_umis": [1], "dna_reads": [10]})
    with pytest.raises(ValueError, match="Malformed"):
        obj(pdf, "mpra_umiwise").validate()


def test_coarse_reporter_tsv_preserves_identifiers_and_deduplicates(tmp_path):
    path = tmp_path / "reporter.tsv"
    path.write_text("rep_id\tcell_bc\tcre_id\tumi\n01\t002\tNA\tu1\n01\t002\tNA\tu2\n")
    data = core.scMPRA_data()
    data.set_coarse_reporter(path)
    assert data._coarse_reporter.to_dict("records") == [{"rep_id": "01", "cell_bc": "002", "cre_id": "NA"}]


@pytest.mark.parametrize("extra", [
    {"transfection_umis": [0]}, {"transfection_reads": [None]},
    {"umis_transfection_bc": [-1]}, {"cell_bc": [None]}, {"cre_id": [""]},
])
def test_coarse_reporter_rejects_invalid_detections(extra):
    reporter = pd.DataFrame({"rep_id": ["r"], "cell_bc": ["c"], "cre_id": ["e"], **extra})
    with pytest.raises(AssertionError):
        core.scMPRA_data().set_coarse_reporter(reporter)


@pytest.mark.parametrize("old_names", [False, True])
@pytest.mark.parametrize("has_reporter", [False, True])
def test_simulation_worker_reads_both_schemas(tmp_path, old_names, has_reporter):
    pdf = pd.DataFrame({"cell_type": ["ct"] * 8,
                        "cre_id": ["reference"] * 4 + ["active"] * 4,
                        "mpra_umis": [0, 1, 2, 1, 0, 8, 9, 10]})
    source = tmp_path / "sim.scmpra"
    source.mkdir()
    (legacy(pdf) if old_names else pdf).to_parquet(source / "data.parquet", index=False)
    hs = core.HypothesisSet.from_dataframe(pd.DataFrame({
        "comparison_CRE": ["active"], "comparison_cell_type": ["ct"],
        "reference_CRE": ["reference"], "reference_cell_type": ["ct"],
    }))
    dest = tmp_path / "result.tsv"
    core._counts_test_worker(None, source, dest, hs, has_reporter,
                             core._mwu_row_fn, "mwu")
    actual = pd.read_csv(dest, sep="\t")
    expected_pdf = pdf if has_reporter else pdf[pdf["mpra_umis"] > 0]
    expected = core._mwu_row_fn(hs.to_dataframe().iloc[0].to_dict(),
                                core._build_mwu_counts_dict(expected_pdf))
    assert actual.loc[0, "p_value"] == pytest.approx(expected["p_value"])
    assert actual.loc[0, "fold_change"] == pytest.approx(expected["fold_change"])


def test_future_schema_version_is_rejected(tmp_path):
    data = obj(readwise())
    dest = tmp_path / "future.scmpra"
    data.to_parquet(dest)
    path = dest / "members.json"
    members = json.loads(path.read_text())
    members["schema_version"] = 999
    path.write_text(json.dumps(members))
    with pytest.raises(ValueError, match="Unsupported MPRA schema_version"):
        core.scMPRA_data.from_parquet(dest)
