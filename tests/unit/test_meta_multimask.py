"""Several masks per gene in one meta run, combined by ACAT (#196)."""

from __future__ import annotations

import math

import numpy as np
import pandas as pd
import pytest

from tests.unit.meta_helpers import make_stratum
from variantcentrifuge.association.correction import apply_correction
from variantcentrifuge.association.meta.artifact import save_artifact
from variantcentrifuge.association.meta.cli import main
from variantcentrifuge.association.meta.manifest import (
    MetaCompatibilityError,
    read_study_manifest,
)
from variantcentrifuge.association.meta.runner import OUTPUT_COLUMNS, run_meta

KEYS = {
    "SIG": [f"1:{100 + j}:A:G" for j in range(12)],
    "NULLG": [f"2:{100 + j}:A:G" for j in range(8)],
    "THIRD": [f"3:{100 + j}:A:G" for j in range(6)],
}
# The "missense" mask is a different variant subset (and so a different mask hash).
MIS_KEYS = {g: keys[: len(keys) // 2] for g, keys in KEYS.items()}


def _strata(mask, keys, seeds=(1, 2)):
    return [
        make_stratum(
            sid,
            seed=seed,
            gene_keys=keys,
            effect_gene="SIG",
            effect_sign=1.0,
            mask_id=mask,
            mask_hash=f"hash-{mask}",
        )
        for sid, seed in zip(("A", "B"), seeds, strict=True)
    ]


def _write(tmp_path, rows, name="study.tsv", with_mask=True):
    header = "stratum_id\tartifact_path" + ("\tmask_id" if with_mask else "")
    lines = [header]
    for mask, art in rows:
        rel = f"{mask}_{art.stratum_id}"
        save_artifact(art, str(tmp_path / rel))
        lines.append(f"{art.stratum_id}\t{rel}" + (f"\t{mask}" if with_mask else ""))
    path = tmp_path / name
    path.write_text("\n".join(lines) + "\n")
    return str(path)


def _two_mask_study(tmp_path):
    rows = [("pLoF", a) for a in _strata("pLoF", KEYS)]
    rows += [("missense", a) for a in _strata("missense", MIS_KEYS, seeds=(3, 4))]
    return rows


def _closed_form_acat(ps):
    """Independent Cauchy combination (Liu & Xie 2020), equal weights, closed form."""
    t = sum(math.tan((0.5 - p) * math.pi) for p in ps) / len(ps)
    return 0.5 - math.atan(t) / math.pi


def test_manifest_mask_column_groups_rows(tmp_path):
    path = _write(tmp_path, _two_mask_study(tmp_path))
    refs = read_study_manifest(path)
    assert [(r.mask_id, r.stratum_id) for r in refs] == [
        ("pLoF", "A"),
        ("pLoF", "B"),
        ("missense", "A"),
        ("missense", "B"),
    ]


def test_manifest_mask_needs_two_strata_and_unique_ids(tmp_path):
    a, b = _strata("pLoF", KEYS)
    path = _write(tmp_path, [("pLoF", a), ("pLoF", b), ("mis", a)])
    with pytest.raises(MetaCompatibilityError, match="mis"):
        read_study_manifest(path)
    path = _write(tmp_path, [("pLoF", a), ("pLoF", a), ("pLoF", b)], name="dup.tsv")
    with pytest.raises(MetaCompatibilityError, match="Duplicate"):
        read_study_manifest(path)


def test_each_mask_matches_its_single_mask_run(tmp_path):
    """Burden/SKAT/SKAT-O per mask are identical to running that mask alone."""
    rows = _two_mask_study(tmp_path)
    multi = run_meta(_write(tmp_path, rows), str(tmp_path / "multi.tsv"))
    assert list(multi.columns) == [*OUTPUT_COLUMNS, "acat_masks_p", "acat_masks_q"]
    for mask in ("pLoF", "missense"):
        single = run_meta(
            _write(tmp_path, [r for r in rows if r[0] == mask], name=f"{mask}.tsv"),
            str(tmp_path / f"{mask}.out"),
        )
        got = multi[multi["mask"] == mask].reset_index(drop=True)
        pd.testing.assert_frame_equal(got[OUTPUT_COLUMNS], single, check_dtype=False)


def test_acat_masks_p_is_cauchy_combination_and_corrected(tmp_path):
    df = run_meta(_write(tmp_path, _two_mask_study(tmp_path)), str(tmp_path / "o.tsv"))
    per_gene = df.groupby("gene")
    acat = per_gene["acat_masks_p"].first()
    for gene, rows in per_gene:
        expected = _closed_form_acat(rows["primary_p"].tolist())
        assert acat[gene] == pytest.approx(expected, rel=1e-12)
        assert rows["acat_masks_p"].nunique() == 1  # repeated on every mask row
    # Family: one ACAT p per gene, corrected across genes with the run's method.
    q = apply_correction(acat.to_numpy(dtype=float), "fdr")
    np.testing.assert_allclose(per_gene["acat_masks_q"].first().to_numpy(), q, rtol=1e-12)


def test_acat_known_value():
    """Hand-checkable value: p=(0.01, 0.2) -> T = (tan(.49π)+tan(.3π))/2."""
    from variantcentrifuge.association.tests.acat import cauchy_combination

    # T = 16.598448937 -> p = 0.5 - atan(T)/π = 0.0191539346 (Python math, closed form)
    assert cauchy_combination([0.01, 0.2]) == pytest.approx(0.0191539346, abs=1e-10)


def test_acat_skips_masks_outside_the_correction_family(tmp_path):
    rows = _two_mask_study(tmp_path)
    df = run_meta(_write(tmp_path, rows), str(tmp_path / "o.tsv"), min_mac=1e9)
    assert df["acat_masks_p"].isna().all() and df["acat_masks_q"].isna().all()


def test_single_mask_output_unchanged_and_label_names_mask(tmp_path):
    a, b = _strata("m", KEYS)
    a.mask_label = b.mask_label = "pLoF"
    df = run_meta(_write(tmp_path, [("x", a), ("x", b)], with_mask=False), str(tmp_path / "o"))
    assert list(df.columns) == OUTPUT_COLUMNS
    assert set(df["mask"]) == {"pLoF"}


def test_mask_groups_validated_separately(tmp_path):
    a, b = _strata("pLoF", KEYS)
    b.mask_hash = "other"
    rows = [("pLoF", a), ("pLoF", b)] + [("mis", x) for x in _strata("mis", MIS_KEYS)]
    with pytest.raises(MetaCompatibilityError, match="mask_hash"):
        run_meta(_write(tmp_path, rows), str(tmp_path / "o"))


def test_cli_multi_mask_diagnostics_per_mask(tmp_path):
    path = _write(tmp_path, _two_mask_study(tmp_path))
    diag = tmp_path / "diag"
    out = tmp_path / "o.tsv"
    assert main(["--manifest", path, "--output", str(out), "--diagnostics-output", str(diag)]) == 0
    for mask in ("pLoF", "missense"):
        lam = pd.read_csv(diag / mask / "lambda_gc.tsv", sep="\t")
        assert set(lam["n_tests"]) == {3}
    assert "acat_masks_p" in pd.read_csv(out, sep="\t").columns


def test_manifest_empty_mask_id_refused(tmp_path):
    m = tmp_path / "m.tsv"
    m.write_text("stratum_id\tartifact_path\tmask_id\nA\ta\tpLoF\nB\tb\t\n")
    with pytest.raises(MetaCompatibilityError, match="empty fields"):
        read_study_manifest(str(m))


def test_acat_masks_keeps_null_mask(tmp_path):
    """A mask with primary_p == 1 pulls acat_masks_p to 1 instead of being dropped."""
    from variantcentrifuge.association.meta.runner import add_acat_masks

    df = pd.DataFrame({"gene": ["G", "G"], "primary_p": [0.01, 1.0], "primary_q": [0.02, 1.0]})
    out = add_acat_masks(df, "fdr")
    assert out["acat_masks_p"].iloc[0] == pytest.approx(1.0, abs=1e-12)


@pytest.mark.parametrize("bad", ["../x", "a/b", ".", "..", "p LoF"])
def test_manifest_mask_id_must_be_a_safe_name(tmp_path, bad):
    """mask_id names a diagnostics subdirectory, so it must not escape DIR."""
    m = tmp_path / "m.tsv"
    m.write_text(f"stratum_id\tartifact_path\tmask_id\nA\ta\t{bad}\nB\tb\t{bad}\n")
    with pytest.raises(MetaCompatibilityError, match="mask_id"):
        read_study_manifest(str(m))
