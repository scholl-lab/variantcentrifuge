"""Unit tests for workflow/rules/helpers.py.

These tests import the helpers module directly — Snakemake is **not** required.
"""

from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

# Make the helpers module importable without Snakemake installed.
_HELPERS_DIR = Path(__file__).resolve().parents[2] / "workflow" / "rules"
if str(_HELPERS_DIR) not in sys.path:
    sys.path.insert(0, str(_HELPERS_DIR))

from helpers import build_vc_flags, get_samples, get_vcf_path


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------
@pytest.fixture()
def samples_df() -> pd.DataFrame:
    """Minimal sample sheet DataFrame."""
    return pd.DataFrame(
        {
            "sample": ["bravo", "alpha", "charlie", "alpha"],
            "vcf_basename": ["bravo_v1", "alpha_v2", "charlie_v3", "alpha_v2"],
        }
    ).set_index("sample", drop=False)


@pytest.fixture()
def minimal_vc_config() -> dict:
    """Minimal variantcentrifuge config dict (all optional flags off)."""
    return {
        "config_file": "config.json",
        "genes": "all",
        "log_level": "INFO",
        "threads": 4,
        "add_chr": False,
        "generate_xlsx": False,
        "generate_html_report": False,
        "field_profile": None,
        "presets": [],
        "append_genotype_fields": [],
        "gene_list_files": [],
    }


@pytest.fixture()
def full_vc_config() -> dict:
    """Fully-populated variantcentrifuge config dict."""
    return {
        "config_file": "config.json",
        "genes": "all",
        "log_level": "INFO",
        "threads": 8,
        "add_chr": True,
        "generate_xlsx": True,
        "generate_html_report": True,
        "field_profile": "dbnsfp5",
        "presets": ["rare", "coding"],
        "append_genotype_fields": ["GEN[*].DP", "GEN[*].AD"],
        "gene_list_files": ["/lists/cardiac.txt", "/lists/cancer.txt"],
    }


@pytest.fixture()
def igv_config() -> dict:
    """IGV config with all options set."""
    return {
        "enabled": True,
        "reference": "hg38_1kg",
        "bam_mapping_file": "/bams/mapping.tsv",
        "fasta": "/ref/genome.fa",
        "ideogram": "/ref/ideogram.json",
        "flanking": 100,
    }


# ---------------------------------------------------------------------------
# get_samples
# ---------------------------------------------------------------------------
class TestGetSamples:
    def test_returns_sorted_unique(self, samples_df: pd.DataFrame) -> None:
        result = get_samples(samples_df)
        assert result == ["alpha", "bravo", "charlie"]

    def test_empty_df(self) -> None:
        df = pd.DataFrame({"sample": [], "vcf_basename": []}).set_index("sample", drop=False)
        assert get_samples(df) == []


# ---------------------------------------------------------------------------
# get_vcf_path
# ---------------------------------------------------------------------------
class TestGetVcfPath:
    def test_resolves_correct_path(self, samples_df: pd.DataFrame) -> None:
        path = get_vcf_path("bravo", "/data/vcfs", samples_df)
        assert path.replace("\\", "/") == "/data/vcfs/bravo_v1.vcf.gz"

    def test_missing_sample_raises(self, samples_df: pd.DataFrame) -> None:
        with pytest.raises(KeyError):
            get_vcf_path("nonexistent", "/data/vcfs", samples_df)


# ---------------------------------------------------------------------------
# build_vc_flags
# ---------------------------------------------------------------------------
class TestBuildVcFlags:
    def test_minimal_config_no_flags(self, minimal_vc_config: dict) -> None:
        flags = build_vc_flags(minimal_vc_config, {})
        assert flags == []

    def test_full_config_all_flags(self, full_vc_config: dict) -> None:
        flags = build_vc_flags(full_vc_config, {})
        assert "--add-chr" in flags
        assert "--xlsx" in flags
        assert "--html-report" in flags
        assert "--field-profile dbnsfp5" in flags
        assert "--preset rare" in flags
        assert "--preset coding" in flags
        assert "--append-extra-sample-fields GEN[*].DP GEN[*].AD" in flags
        assert "--annotate-gene-list /lists/cardiac.txt" in flags
        assert "--annotate-gene-list /lists/cancer.txt" in flags

    def test_igv_flags(self, minimal_vc_config: dict, igv_config: dict) -> None:
        flags = build_vc_flags(minimal_vc_config, igv_config)
        assert "--igv" in flags
        assert "--igv-reference hg38_1kg" in flags
        assert "--bam-mapping-file /bams/mapping.tsv" in flags
        assert "--igv-fasta /ref/genome.fa" in flags
        assert "--igv-ideogram /ref/ideogram.json" in flags
        assert "--igv-flanking 100" in flags

    def test_igv_disabled_no_igv_flags(self, minimal_vc_config: dict) -> None:
        igv = {"enabled": False, "reference": "hg19"}
        flags = build_vc_flags(minimal_vc_config, igv)
        assert not any("--igv" in f for f in flags)

    def test_field_profile_none_skipped(self, minimal_vc_config: dict) -> None:
        flags = build_vc_flags(minimal_vc_config, {})
        assert not any("--field-profile" in f for f in flags)

    def test_igv_default_flanking(self, minimal_vc_config: dict) -> None:
        igv = {"enabled": True}
        flags = build_vc_flags(minimal_vc_config, igv)
        assert "--igv-flanking 50" in flags
