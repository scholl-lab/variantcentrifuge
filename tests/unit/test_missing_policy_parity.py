"""Genotype missing-data policies vs R SKAT 2.2.5 (#200).

Fixtures and references come from ``scripts/generate_missing_impute_parity.R``
(R 4.5.3, SKAT 2.2.5; provenance line at the top of ``r_reference.csv``). R read the
same CSVs (NA = missing call) and ran ``SKAT(..., method="davies", missing_cutoff=0.15)``
with ``impute.method="bestguess"`` (binary; the SKATBinary default) or ``"fixed"``
(quantitative; the SKAT default), and with ``Z[is.na(Z)] <- 0`` for ``missing_as_ref``.
One variant (g9) is 30% missing: both implementations drop it under the default policy
and keep it (as REF calls) under ``missing_as_ref``.

Tolerance ``rel=1e-6`` as in ``test_skat_null_parity.py`` (Davies ``acc=1e-6``).
"""

from __future__ import annotations

import csv
import pathlib

import numpy as np
import pandas as pd
import pytest

from variantcentrifuge.association.backends.python_backend import PythonSKATBackend
from variantcentrifuge.association.genotype_matrix import build_genotype_matrix

FIXTURES = pathlib.Path(__file__).parent.parent / "fixtures" / "missing_impute_parity"
GT = {0.0: "0/0", 1.0: "0/1", 2.0: "1/1"}


def _references() -> dict[tuple[str, str], float]:
    with open(FIXTURES / "r_reference.csv") as fh:
        rows = csv.DictReader(line for line in fh if not line.startswith("#"))
        return {(r["name"], r["method"]): float(r["p_value"]) for r in rows}


R_REF = _references()


@pytest.fixture(scope="module")
def backend() -> PythonSKATBackend:
    b = PythonSKATBackend()
    b.detect_environment()
    return b


@pytest.mark.parametrize("missing_as_ref", [False, True], ids=["default", "missing_as_ref"])
@pytest.mark.parametrize("method", ["SKAT", "Burden"])
@pytest.mark.parametrize("trait", ["binary", "quantitative"])
def test_missing_policy_matches_r(backend, trait, method, missing_as_ref):
    data = pd.read_csv(FIXTURES / f"{trait}.csv")
    dosage = data.filter(regex=r"^g\d+$").to_numpy(dtype=float)  # (n_samples, p), NaN = missing
    samples = [f"S{i}" for i in range(dosage.shape[0])]
    gene_df = pd.DataFrame(
        {
            s: ["./." if np.isnan(d) else GT[d] for d in row]
            for s, row in zip(samples, dosage, strict=True)
        }
    )
    with open(FIXTURES / "weights.csv") as fh:
        w = np.array([float(r["w"]) for r in csv.DictReader(fh) if r["name"] == trait])

    geno, _mafs, mask, _warn, keep = build_genotype_matrix(  # type: ignore[misc]
        gene_df,
        samples,
        samples,
        is_binary=trait == "binary",
        missing_site_threshold=0.15,
        missing_as_ref=missing_as_ref,
        return_keep_mask=True,
    )
    assert all(mask)
    assert bool(keep[8]) is missing_as_ref  # g9 (30% missing)
    null = backend.fit_null_model(data["y"].to_numpy(dtype=float), None, trait)
    p = backend.test_gene("G", geno, null, method, weights=w[keep])["p_value"]
    policy = "missing_as_ref" if missing_as_ref else "default"
    assert p == pytest.approx(R_REF[(trait, f"{method}_{policy}")], rel=1e-6, abs=0)
