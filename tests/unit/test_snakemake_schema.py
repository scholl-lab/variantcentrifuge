"""Unit tests for the Snakemake workflow schema files.

Validates the JSON Schema documents under ``workflow/schemas/`` against
sample config and sample-sheet data without requiring Snakemake.
"""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

try:
    import jsonschema
except ImportError:  # pragma: no cover
    jsonschema = None  # type: ignore[assignment]

SCHEMA_DIR = Path(__file__).resolve().parents[2] / "workflow" / "schemas"

skip_if_no_jsonschema = pytest.mark.skipif(
    jsonschema is None,
    reason="jsonschema package not installed",
)


def _load_schema(name: str) -> dict:
    with open(SCHEMA_DIR / name, encoding="utf-8") as fh:
        return yaml.safe_load(fh)


# ---------------------------------------------------------------------------
# Config schema
# ---------------------------------------------------------------------------
@skip_if_no_jsonschema
class TestConfigSchema:
    @pytest.fixture()
    def schema(self) -> dict:
        return _load_schema("config.schema.yaml")

    @pytest.fixture()
    def valid_config(self) -> dict:
        return {
            "paths": {
                "samples": "config/samples.tsv",
                "vcf_folder": "/data/vcfs",
                "output_folder": "results/vc",
            },
            "variantcentrifuge": {
                "config_file": "/etc/vc/config.json",
            },
        }

    def test_valid_minimal_config(self, schema: dict, valid_config: dict) -> None:
        jsonschema.validate(valid_config, schema)

    def test_missing_paths_vcf_folder(self, schema: dict) -> None:
        bad = {
            "paths": {
                "samples": "s.tsv",
                "output_folder": "out",
                # vcf_folder missing
            },
            "variantcentrifuge": {"config_file": "c.json"},
        }
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(bad, schema)

    def test_missing_config_file(self, schema: dict) -> None:
        bad = {
            "paths": {
                "samples": "s.tsv",
                "vcf_folder": "/v",
                "output_folder": "out",
            },
            "variantcentrifuge": {},
        }
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(bad, schema)

    def test_invalid_log_level(self, schema: dict, valid_config: dict) -> None:
        valid_config["variantcentrifuge"]["log_level"] = "VERBOSE"
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate(valid_config, schema)

    def test_full_config(self, schema: dict) -> None:
        full = {
            "paths": {
                "samples": "config/samples.tsv",
                "vcf_folder": "/data/vcfs",
                "output_folder": "results/vc",
                "log_subdir": "logs",
            },
            "variantcentrifuge": {
                "config_file": "config.json",
                "genes": "all",
                "log_level": "DEBUG",
                "threads": 4,
                "add_chr": True,
                "generate_xlsx": True,
                "generate_html_report": False,
                "field_profile": "dbnsfp5",
                "presets": ["rare", "coding"],
                "append_genotype_fields": ["GEN[*].DP"],
                "gene_list_files": ["/lists/cancer.txt"],
            },
            "igv": {
                "enabled": True,
                "reference": "hg19",
                "bam_mapping_file": "bams.tsv",
                "fasta": "/ref/genome.fa",
                "ideogram": None,
                "flanking": 100,
            },
            "container": {
                "image": "docker://ghcr.io/scholl-lab/variantcentrifuge:0.10.0",
            },
        }
        jsonschema.validate(full, schema)


# ---------------------------------------------------------------------------
# Samples schema
# ---------------------------------------------------------------------------
@skip_if_no_jsonschema
class TestSamplesSchema:
    @pytest.fixture()
    def schema(self) -> dict:
        return _load_schema("samples.schema.yaml")

    def test_valid_sample(self, schema: dict) -> None:
        jsonschema.validate(
            {"sample": "sample_01", "vcf_basename": "sample_01_filtered"},
            schema,
        )

    def test_missing_sample_field(self, schema: dict) -> None:
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate({"vcf_basename": "test"}, schema)

    def test_missing_vcf_basename(self, schema: dict) -> None:
        with pytest.raises(jsonschema.ValidationError):
            jsonschema.validate({"sample": "test"}, schema)
