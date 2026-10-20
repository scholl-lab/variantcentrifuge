"""Real-data e2e: two 1000 Genomes cohorts -> --score-artifact -> variantcentrifuge meta.

Fixture: tests/fixtures/1000g_cohort (scripts/build_1000g_fixture.sh). Real EUR and SAS
genotypes for 10 nephrology genes, annotated like AGDE (snpEff 5.2 GRCh38.p14 -lof,
ClinVar 2025-05-04); carriers of rare HIGH/MODERATE PKHD1 variants are seeded as likely
cases, so PKHD1 is the one true signal. Needs bcftools and SnpSift: `make test-e2e`.
"""

from __future__ import annotations

import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

from variantcentrifuge.association.meta.artifact import load_artifact

FIXTURE = Path(__file__).parents[1] / "fixtures" / "1000g_cohort"
VCF = FIXTURE / "1000g_nephro.GRCh38.annotated.vcf.gz"

pytestmark = [
    pytest.mark.e2e,
    pytest.mark.slow,
    pytest.mark.skipif(
        not (shutil.which("bcftools") and shutil.which("SnpSift")),
        reason="needs bcftools and SnpSift (run: make test-e2e)",
    ),
]


def _vc(*args: str) -> None:
    cmd = [
        sys.executable,
        "-c",
        "import sys; from variantcentrifuge.cli import main; sys.exit(main())",
    ]
    res = subprocess.run([*cmd, *args], capture_output=True, text=True, timeout=600)
    assert res.returncode == 0, res.stderr[-3000:]


@pytest.fixture(scope="module")
def study(tmp_path_factory: pytest.TempPathFactory) -> Path:
    out = tmp_path_factory.mktemp("1000g")
    for cohort in ("EUR", "SAS"):
        _vc(
            "-v", str(VCF),
            "--reference", "GRCh38.p14",
            "--gene-name", "all",
            "--regions-bed", str(FIXTURE / "genes.GRCh38.bed"),
            "--split-snpeff-lines", "before_filters",
            "--bcftools-prefilter", "INFO/AF < 0.01",
            "--filters", "((ANN[0].IMPACT has 'HIGH') | (ANN[0].IMPACT has 'MODERATE'))",
            "--fields", "CHROM POS REF ALT AF ANN[0].GENE ANN[0].EFFECT ANN[0].IMPACT "
            "ANN[0].FEATUREID CLNSIG GEN[*].GT",
            "--case-samples-file", str(FIXTURE / f"{cohort}.cases.txt"),
            "--control-samples-file", str(FIXTURE / f"{cohort}.controls.txt"),
            "--perform-association",
            "--association-tests", "fisher,logistic_burden,skat",
            "--score-artifact", str(out / f"{cohort}.npz"),
            "--output-dir", str(out / cohort),
            "--output-file", f"{cohort}.tsv",
            "--no-stats",
        )  # fmt: skip
    (out / "study.tsv").write_text("stratum_id\tartifact_path\nEUR\tEUR.npz\nSAS\tSAS.npz\n")
    _vc("meta", "--manifest", str(out / "study.tsv"), "--output", str(out / "meta.tsv"))
    return out


def test_artifacts_match_cohorts(study: Path) -> None:
    arts = {c: load_artifact(str(study / f"{c}.npz")) for c in ("EUR", "SAS")}
    for cohort, art in arts.items():
        n_cases = len((FIXTURE / f"{cohort}.cases.txt").read_text().split())
        n_controls = len((FIXTURE / f"{cohort}.controls.txt").read_text().split())
        assert (art.n_cases, art.n_controls) == (n_cases, n_controls)
        assert art.genome_build == "GRCh38" and "PKHD1" in art.genes
    assert not set(arts["EUR"].sample_id_hashes) & set(arts["SAS"].sample_id_hashes)


def test_meta_recovers_planted_gene_only(study: Path) -> None:
    meta = pd.read_csv(study / "meta.tsv", sep="\t").set_index("gene")
    pkhd1 = meta.loc["PKHD1"]
    assert pkhd1["n_strata"] == 2 and pkhd1["status"] == "ok"
    assert pkhd1["primary_q"] < 1e-10
    assert pkhd1["per_stratum_directions"] == "EUR:1,SAS:1"
    others = meta.drop(index="PKHD1")
    assert (others["primary_q"] > 0.05).all(), others[["primary_p", "primary_q"]]
    # Every gene got a result: one-variant genes fall back to burden instead of crashing.
    assert meta["primary_p"].notna().all()
