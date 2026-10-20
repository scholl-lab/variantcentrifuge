"""Seeded synthetic case/control labels for the 1000 Genomes fixture (see build_1000g_fixture.sh).

Carriers of a rare (AF < 1%) HIGH or MODERATE variant in the signal gene are cases with
probability 0.6, everyone else with 0.2, so the signal gene is the one true association and
the other genes are real-genotype nulls. Writes <out>/<cohort>.cases.txt / .controls.txt.
"""

from __future__ import annotations

import gzip
import sys
from pathlib import Path

import numpy as np

SIGNAL_GENE = "PKHD1"
P_CASE_CARRIER, P_CASE_OTHER = 0.6, 0.2
SEED = 112  # issue #112


def carriers(vcf: Path) -> set[str]:
    found: set[str] = set()
    samples: list[str] = []
    with gzip.open(vcf, "rt") as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            cols = line.rstrip("\n").split("\t")
            if line.startswith("#"):
                samples = cols[9:]
                continue
            info = dict(kv.split("=", 1) for kv in cols[7].split(";") if "=" in kv)
            if float(info["AF"]) >= 0.01:
                continue
            ann = [a.split("|") for a in info.get("ANN", "").split(",")]
            if not any(a[3] == SIGNAL_GENE and a[2] in ("HIGH", "MODERATE") for a in ann):
                continue
            found.update(s for s, gt in zip(samples, cols[9:], strict=True) if "1" in gt)
    return found


def main(out: Path) -> None:
    hits = carriers(out / "1000g_nephro.GRCh38.annotated.vcf.gz")
    rows = [line.split("\t") for line in (out / "populations.tsv").read_text().splitlines()]
    rng = np.random.default_rng(SEED)
    for cohort in ("EUR", "SAS"):
        ids = [s for s, pop in rows if pop == cohort]
        p = np.array([P_CASE_CARRIER if s in hits else P_CASE_OTHER for s in ids])
        is_case = rng.uniform(size=len(ids)) < p
        cases = [s for s, c in zip(ids, is_case, strict=True) if c]
        controls = [s for s, c in zip(ids, is_case, strict=True) if not c]
        (out / f"{cohort}.cases.txt").write_text("\n".join(cases) + "\n")
        (out / f"{cohort}.controls.txt").write_text("\n".join(controls) + "\n")
        n_carriers = sum(s in hits for s in ids)
        print(f"{cohort}: {len(cases)} cases, {len(controls)} controls, {n_carriers} carriers")


if __name__ == "__main__":
    main(Path(sys.argv[1]))
