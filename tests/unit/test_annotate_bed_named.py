"""Named/gz BED annotation and 0-based vs 1-based edge handling (#94)."""

import gzip

import pandas as pd
import pytest

from variantcentrifuge.annotator import annotate_dataframe_with_features, load_custom_features
from variantcentrifuge.validators import validate_annotate_bed_specs

# BED chr1 100 110 = 1-based POS 101..110
BED = "chr1\t100\t110\tA\n"


def _run(specs, rows):
    df = pd.DataFrame(rows, columns=["CHROM", "POS", "REF", "ALT"])
    feats = load_custom_features({"annotate_bed_files": specs})
    return annotate_dataframe_with_features(df, feats)


@pytest.fixture
def bed(tmp_path):
    p = tmp_path / "t.bed"
    p.write_text(BED)
    return str(p)


def test_named_edges_and_deletion(bed):
    rows = [
        ("1", 100, "A", "T"),  # before start
        ("1", 101, "A", "T"),  # first base
        ("1", 110, "A", "T"),  # last base
        ("1", 111, "A", "T"),  # after end
        ("1", 95, "A" * 12, "A"),  # deletion spanning into interval
    ]
    out = _run([f"Rep={bed}"], rows)
    assert list(out["Rep"]) == ["no", "yes", "yes", "no", "yes"]
    assert (out["Custom_Annotation"] == "").all()


def test_legacy_edges_and_gz_with_named(tmp_path, bed):
    gz = tmp_path / "g.bed.gz"
    with gzip.open(gz, "wt") as f:
        f.write(BED)
    out = _run([bed, f"Gz={gz}"], [("chr1", 100, "A", "T"), ("chr1", 110, "A", "T")])
    assert list(out["Custom_Annotation"]) == ["", "Region=A_t"]
    assert list(out["Gz"]) == ["no", "yes"]


def test_unreadable_raises(tmp_path):
    bad = tmp_path / "bad.bed.gz"
    bad.write_text("not gzip")
    with pytest.raises(OSError):
        load_custom_features({"annotate_bed_files": [f"X={bad}"]})


def test_column_clash_raises(bed):
    df = pd.DataFrame({"CHROM": ["1"], "POS": [101], "REF": ["A"], "Rep": ["z"]})
    feats = load_custom_features({"annotate_bed_files": [f"Rep={bed}"]})
    with pytest.raises(ValueError, match="clashes"):
        annotate_dataframe_with_features(df, feats)


def test_validate_specs(bed):
    validate_annotate_bed_specs([f"A={bed}", bed])
    with pytest.raises(ValueError):
        validate_annotate_bed_specs([f"A={bed}", f"A={bed}"])
    with pytest.raises(ValueError):
        validate_annotate_bed_specs([f"GT={bed}"])
    with pytest.raises(FileNotFoundError):
        validate_annotate_bed_specs(["A=/no/such.bed"])


def test_checkpoint_inputs_track_named_bed(bed):
    from types import SimpleNamespace

    from variantcentrifuge.stages.data_loading_stages import CustomAnnotationStage

    ctx = SimpleNamespace(annotation_configs={"bed_files": [f"Rep={bed}"]})
    assert [str(p) for p in CustomAnnotationStage().get_input_files(ctx)] == [bed]


def test_fields_clash_rejected(bed):
    from types import SimpleNamespace

    from variantcentrifuge.validators import validate_cli_args

    args = SimpleNamespace(annotate_bed=[f"IMPACT={bed}"], fields="CHROM POS ANN[0].IMPACT")
    with pytest.raises(ValueError, match="IMPACT"):
        validate_cli_args(args, {})


def test_nan_pos_and_ref_and_special_alleles(bed):
    rows = [
        ("1", float("nan"), "A", "T"),  # NaN POS -> no
        ("1", 105, float("nan"), "T"),  # NaN REF -> length 1
        ("1", 105, "A", "<DEL>"),  # symbolic ALT: REF span only
        ("1", 105, "*", "A"),
        ("1", 111, ".", "A"),
        ("1", 110, "A", "AT"),  # insertion anchored on last base
        ("1", 100, "A", "AT"),  # insertion anchored before start
        ("MT", 105, "A", "T"),  # contig absent from BED
    ]
    out = _run([f"Rep={bed}"], rows)
    assert list(out["Rep"]) == ["no", "yes", "yes", "yes", "no", "yes", "no", "no"]


def test_chr_mismatch_matches_legacy(tmp_path):
    p = tmp_path / "m.bed"
    p.write_text("chrM\t100\t110\tA\n")
    out = _run([str(p), f"N={p}"], [("MT", 105, "A", "T"), ("chrM", 105, "A", "T")])
    assert list(out["N"]) == ["no", "yes"]
    assert list(out["Custom_Annotation"]) == ["", "Region=A_m"]
