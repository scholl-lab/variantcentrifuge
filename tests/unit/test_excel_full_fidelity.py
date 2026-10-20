"""Excel output fidelity test: internal cache columns are excluded from output."""

import numpy as np
import pandas as pd
import pytest
from openpyxl import load_workbook

from variantcentrifuge.converter import convert_to_excel


def _create_test_variants_df(n_variants: int = 20, n_samples: int = 3, seed: int = 42):
    """Generate synthetic variant DataFrame for testing."""
    rng = np.random.default_rng(seed)

    # Generate genomic data
    chroms = rng.choice(["chr1", "chr2", "chr3"], size=n_variants)
    positions = rng.integers(1000, 100000, size=n_variants)
    refs = rng.choice(["A", "C", "G", "T"], size=n_variants)
    alts = rng.choice(["A", "C", "G", "T"], size=n_variants)

    # Generate GT column
    gts = []
    for _ in range(n_variants):
        sample_gts = [
            f"Sample{i + 1}({rng.choice(['0/0', '0/1', '1/1', './.'])})" for i in range(n_samples)
        ]
        gts.append(";".join(sample_gts))

    # Generate annotations and URL columns
    genes = rng.choice(["BRCA1", "TP53", "EGFR"], size=n_variants)
    impacts = rng.choice(["HIGH", "MODERATE", "LOW"], size=n_variants)

    # Add URL columns for hyperlink testing
    splicai_urls = [
        f"https://spliceailookup.broadinstitute.org/#variant={chrom}-{pos}-{ref}-{alt}"
        for chrom, pos, ref, alt in zip(chroms, positions, refs, alts, strict=True)
    ]
    gnomad_urls = [
        f"https://gnomad.broadinstitute.org/variant/{chrom}-{pos}-{ref}-{alt}"
        for chrom, pos, ref, alt in zip(chroms, positions, refs, alts, strict=True)
    ]

    return pd.DataFrame(
        {
            "CHROM": chroms,
            "POS": positions,
            "REF": refs,
            "ALT": alts,
            "GT": gts,
            "GENE": genes,
            "IMPACT": impacts,
            "SpliceAI_URL": splicai_urls,
            "gnomAD_URL": gnomad_urls,
        }
    )


@pytest.mark.unit
def test_gt_cache_cleanup_before_output(tmp_path):
    """
    Verify cache column cleanup pattern works correctly.

    Tests the cache cleanup pattern used in output stages:
    Drop all columns starting with "_" before writing to Excel.

    NOTE: convert_to_excel itself doesn't strip cache columns - that's
    the responsibility of the calling stage (TSVOutputStage, ExcelReportStage).
    This test verifies the cleanup pattern works correctly.
    """
    # Create test data with cache columns
    variants_df = _create_test_variants_df(n_variants=10)

    # Add cache column (simulating internal processing cache)
    # Note: _GT_PARSED removed in Phase 11 (dead code elimination)
    variants_df["_INTERNAL_CACHE"] = ["cache_value"] * len(variants_df)

    # Apply the cache cleanup pattern from output stages
    cache_cols = [c for c in variants_df.columns if c.startswith("_")]
    assert len(cache_cols) == 1, f"Expected 1 cache column, found {len(cache_cols)}"

    df_cleaned = variants_df.drop(columns=cache_cols)

    # Write TSV and generate Excel with cleaned DataFrame
    variants_tsv = tmp_path / "variants.tsv"
    df_cleaned.to_csv(variants_tsv, sep="\t", index=False)

    cfg = {"links": {}, "igv_enabled": False}
    xlsx_path = convert_to_excel(str(variants_tsv), cfg, df=df_cleaned)

    # Verify no cache columns in Excel
    wb = load_workbook(xlsx_path)
    ws = wb["Results"]
    header_row = [cell.value for cell in ws[1]]

    # Verify all expected columns are present
    expected_cols = [
        "CHROM",
        "POS",
        "REF",
        "ALT",
        "GT",
        "GENE",
        "IMPACT",
        "SpliceAI_URL",
        "gnomAD_URL",
    ]
    for col in expected_cols:
        assert col in header_row, f"Expected column '{col}' missing from output"

    # Verify no cache columns (starting with _)
    cache_cols_in_output = [c for c in header_row if c and c.startswith("_")]
    assert len(cache_cols_in_output) == 0, (
        f"Cache columns should not appear in Excel output. Found: {cache_cols_in_output}"
    )
