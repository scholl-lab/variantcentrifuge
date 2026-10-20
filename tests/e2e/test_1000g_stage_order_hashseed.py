"""Real-data e2e: stage order and association output must not depend on PYTHONHASHSEED (#219).

On v0.22 the EUR cohort run ordered association before variant analysis for seed 0 but
after it for seed 4, switching association aggregation paths. Needs bcftools and SnpSift.
"""

from __future__ import annotations

import gzip
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pytest

FIXTURE = Path(__file__).parents[1] / "fixtures" / "1000g_cohort"

pytestmark = [
    pytest.mark.e2e,
    pytest.mark.slow,
    pytest.mark.skipif(
        not (shutil.which("bcftools") and shutil.which("SnpSift")),
        reason="needs bcftools and SnpSift (run: make test-e2e)",
    ),
]

# Planned levels (DEBUG) and association branch; "Executing" lines of parallel-safe stages
# within one level are logged in completion order, so they are not compared.
_ORDER_LINE = re.compile(r"Level \d+: \[.*|Reconstructing GT column for .*|using .* aggregation")


def _run(seed: str, out: Path) -> tuple[list[str], str, str]:
    cmd = [
        sys.executable, "-c",
        "import sys; from variantcentrifuge.cli import main; sys.exit(main())",
        "-v", str(FIXTURE / "1000g_nephro.GRCh38.annotated.vcf.gz"),
        "--reference", "GRCh38.p14",
        "--gene-name", "all",
        "--regions-bed", str(FIXTURE / "genes.GRCh38.bed"),
        "--bcftools-prefilter", "INFO/AF < 0.01",
        "--filters", "((ANN[0].IMPACT has 'HIGH') | (ANN[0].IMPACT has 'MODERATE'))",
        "--fields", "CHROM POS REF ALT AF ANN[0].GENE ANN[0].EFFECT ANN[0].IMPACT CLNSIG GEN[*].GT",
        "--case-samples-file", str(FIXTURE / "EUR.cases.txt"),
        "--control-samples-file", str(FIXTURE / "EUR.controls.txt"),
        "--perform-association",
        "--association-tests", "fisher,logistic_burden",
        "--output-dir", str(out),
        "--output-file", "out.tsv",
        "--no-stats",
        "--log-level", "DEBUG",
    ]  # fmt: skip
    env = {**os.environ, "PYTHONHASHSEED": seed}
    res = subprocess.run(cmd, capture_output=True, text=True, timeout=600, env=env)
    assert res.returncode == 0, res.stderr[-3000:]
    order = [m.group(0) for m in map(_ORDER_LINE.search, res.stderr.splitlines()) if m]
    assoc = gzip.decompress((out / "association_results.association.tsv.gz").read_bytes())
    return order, assoc.decode(), (out / "out.tsv").read_text()


def test_stage_order_and_association_independent_of_hash_seed(tmp_path: Path) -> None:
    runs = {seed: _run(seed, tmp_path / seed) for seed in ("0", "4")}
    (order0, assoc0, tsv0), (order4, assoc4, tsv4) = runs["0"], runs["4"]
    assert order0 == order4 and any("association_analysis" in line for line in order0)
    assert not any("column-based" in line for line in order0)
    assert assoc0 == assoc4 and assoc0.count("\n") > 1
    assert tsv0 == tsv4


def test_pseudonymize_keeps_proband_counts_and_association(tmp_path: Path) -> None:
    """Codex FX-2: pseudonymizing before variant/association analysis zeroed case counts."""
    import pandas as pd

    def run(out: Path, *extra: str) -> tuple[pd.DataFrame, str]:
        cmd = [
            sys.executable, "-c",
            "import sys; from variantcentrifuge.cli import main; sys.exit(main())",
            "-v", str(FIXTURE / "1000g_nephro.GRCh38.annotated.vcf.gz"),
            "--reference", "GRCh38.p14", "--gene-name", "all",
            "--regions-bed", str(FIXTURE / "genes.GRCh38.bed"),
            "--bcftools-prefilter", "INFO/AF < 0.01",
            "--filters", "((ANN[0].IMPACT has 'HIGH') | (ANN[0].IMPACT has 'MODERATE'))",
            "--fields", "CHROM POS REF ALT ANN[0].GENE ANN[0].IMPACT GEN[*].GT",
            "--case-samples-file", str(FIXTURE / "EUR.cases.txt"),
            "--control-samples-file", str(FIXTURE / "EUR.controls.txt"),
            "--perform-association", "--association-tests", "fisher",
            "--output-dir", str(out), "--output-file", "out.tsv", *extra,
        ]  # fmt: skip
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        assert res.returncode == 0, res.stderr[-3000:]
        assoc = gzip.decompress((out / "association_results.association.tsv.gz").read_bytes())
        return pd.read_csv(out / "out.tsv", sep="\t", dtype=str), assoc.decode()

    plain, assoc_plain = run(tmp_path / "plain")
    pseudo, assoc_pseudo = run(
        tmp_path / "pseudo", "--pseudonymize", "--pseudonymize-table", str(tmp_path / "map.tsv")
    )
    counts = [c for c in plain.columns if c.startswith(("proband_", "control_"))]
    assert counts and plain[counts].astype(float).to_numpy().sum() > 0
    pd.testing.assert_frame_equal(plain[counts], pseudo[counts])
    assert assoc_plain == assoc_pseudo and assoc_plain.count("\n") > 1


def test_default_statistics_identical_with_and_without_sample_subset(tmp_path: Path) -> None:
    """FX-6: statistics no longer block subsetting; stats and association are unchanged."""

    def run(out: Path, *extra: str) -> tuple[str, str, str]:
        cmd = [
            sys.executable, "-c",
            "import sys; from variantcentrifuge.cli import main; sys.exit(main())",
            "-v", str(FIXTURE / "1000g_nephro.GRCh38.annotated.vcf.gz"),
            "--reference", "GRCh38.p14", "--gene-name", "all",
            "--regions-bed", str(FIXTURE / "genes.GRCh38.bed"),
            "--bcftools-prefilter", "INFO/AF < 0.01",
            "--filters", "((ANN[0].IMPACT has 'HIGH') | (ANN[0].IMPACT has 'MODERATE'))",
            "--fields", "CHROM POS REF ALT ANN[0].GENE ANN[0].IMPACT GEN[*].GT",
            "--case-samples-file", str(FIXTURE / "EUR.cases.txt"),
            "--control-samples-file", str(FIXTURE / "EUR.controls.txt"),
            "--perform-association", "--association-tests", "fisher,logistic_burden",
            "--output-dir", str(out), "--output-file", "out.tsv",
            "--stats-output-file", str(out / "stats.tsv"), *extra,
        ]  # fmt: skip
        res = subprocess.run(cmd, capture_output=True, text=True, timeout=600)
        assert res.returncode == 0, res.stderr[-3000:]
        assoc = gzip.decompress((out / "association_results.association.tsv.gz").read_bytes())
        return (out / "stats.tsv").read_text(), assoc.decode(), res.stderr

    stats_sub, assoc_sub, log_sub = run(tmp_path / "subset")
    stats_all, assoc_all, _ = run(tmp_path / "all", "--no-sample-subset")
    assert "Sample subsetting at extraction:" in log_sub
    assert stats_sub == stats_all and "total_samples\t1039" in stats_sub
    assert assoc_sub == assoc_all and assoc_sub.count("\n") > 1
