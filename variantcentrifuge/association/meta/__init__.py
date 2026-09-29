# File: variantcentrifuge/association/meta/__init__.py
"""
Cohort-stratified rare-variant score-statistic meta-analysis (issue #112).

Combines per-stratum score/covariance artifacts (burden, SKAT, SKAT-O) computed
by the existing single-cohort association engine into a common-effect meta-analysis
with heterogeneity diagnostics. Binary traits only (quantitative deferred; see the
design spec revision 1). Reuses the Davies chain and SKAT-O omnibus in
``association.backends`` — nothing is reimplemented that already exists.

Public surface
--------------
- ``ScoreCovarianceArtifact`` / ``GeneScore`` — the persisted aggregate artifact.
- ``compute_gene_score`` / ``fit_stratum_null`` — build score/covariance per gene.
- ``meta_burden`` / ``meta_skat`` / ``meta_skato`` — combine aligned strata.
- ``run_meta`` — end-to-end manifest → meta results TSV.
"""

from __future__ import annotations

from variantcentrifuge.association.meta.artifact import (
    GeneScore,
    ScoreCovarianceArtifact,
    compute_gene_score,
    fit_stratum_null,
    load_artifact,
    save_artifact,
)
from variantcentrifuge.association.meta.combine import meta_burden, meta_skat, meta_skato

__all__ = [
    "GeneScore",
    "ScoreCovarianceArtifact",
    "compute_gene_score",
    "fit_stratum_null",
    "load_artifact",
    "meta_burden",
    "meta_skat",
    "meta_skato",
    "save_artifact",
]
