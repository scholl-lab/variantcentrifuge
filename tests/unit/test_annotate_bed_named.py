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
