"""Restrict field extraction to the analysed samples of an association run (#212).

Like PLINK ``--keep``, REGENIE and SAIGE, an association run with a sample-defining input
(case/control samples or files, a phenotype file) only extracts the per-sample columns of
the analysed set: ``bcftools query -S`` at field extraction. Extraction is the earliest
step that is safe for every filter: upstream bcftools ``-i`` and SnpSift filters still see
all samples, and ``bcftools query`` prints INFO fields as stored (no AC/AN recomputation).
Samples named in a sample file but absent from the VCF are never passed to bcftools, so
they keep today's behaviour (counted in the cohort totals, never genotyped).
"""

from __future__ import annotations

import hashlib
import logging
from argparse import Namespace
from typing import Any

logger = logging.getLogger(__name__)

# Resolved once per run (context attribute, not config: the metadata TSV dumps the config).
_SUBSET_ATTR = "_bcftools_sample_subset"
_VCF_COUNT_ATTR = "_vcf_sample_count_before_subset"


def sample_set_sha256(samples: list[str]) -> str:
    """SHA-256 of the sorted, newline-joined sample IDs."""
    return hashlib.sha256("\n".join(sorted(samples)).encode("utf-8")).hexdigest()


def blocking_reason(cfg: dict[str, Any]) -> str | None:
    """Why all samples must be kept (a feature reads samples outside the analysed set).

    Config only, so the checkpoint key (``hash_configuration``) can include the decision.
    """
    from .scoring_requirements import check_scoring_requires_inheritance

    checks = (
        (
            cfg.get("calculate_inheritance") or cfg.get("ped") or cfg.get("sample_sheet"),
            "inheritance/pedigree analysis",
        ),
        (check_scoring_requires_inheritance(Namespace(), cfg), "scoring needs inheritance"),
        (cfg.get("append_extra_sample_fields"), "--append-extra-sample-fields"),
        (cfg.get("no_replacement"), "--no-replacement (raw per-sample output columns)"),
        (cfg.get("igv_enabled"), "--igv (per-sample IGV reports)"),
        (
            not cfg.get("no_stats") and cfg.get("stats_config"),
            "custom --stats-config (expressions may read per-sample columns)",
        ),
        (cfg.get("pseudonymize"), "--pseudonymize (pseudonyms numbered by sample position)"),
        (
            any(
                cfg.get(k) not in (None, "") for k in ("tumor_sample_index", "normal_sample_index")
            ),
            "tumor/normal sample index (positional sample labels)",
        ),
        (
            cfg.get("late_filtering") or cfg.get("final_filter"),
            "--late-filtering/--final-filter (expressions may read per-sample columns)",
        ),
    )
    return next((reason for active, reason in checks if active), None)


def _analysed_samples(context: Any, vcf_samples: list[str]) -> set[str] | None:
    """The analysed sample set, or None when no input defines one."""
    cfg = context.config
    if cfg.get("trait_type") == "quantitative":
        if not (cfg.get("phenotype_file") and cfg.get("phenotype_value_column")):
            return None
        from .stages.association_outcome import build_quantitative_outcome

        samples, _ = build_quantitative_outcome(
            context.phenotype_data,
            vcf_samples,
            cfg.get("phenotype_value_column"),
            cfg.get("remove_sample_substring"),
        )
        return set(samples)
    cases, controls = cfg.get("case_samples") or [], cfg.get("control_samples") or []
    if isinstance(cases, str) or isinstance(controls, str):
        return None  # unparsed CLI string: the association stage parses and validates it
    return set(cases) | set(controls) if cases and controls else None


def resolve_sample_subset(context: Any) -> list[str] | None:
    """Return raw VCF sample names (VCF order) to extract, or None to keep every sample.

    When a subset applies, ``context.vcf_samples`` and ``config["sample_list"]`` are narrowed
    to it and ``config["sample_subset_sha256"]`` records it. Idempotent: the first decision
    is cached on the context, so extraction, chunks and checkpoint restores agree.
    """
    cfg = context.config
    if _SUBSET_ATTR in vars(context):
        cached: list[str] | None = vars(context)[_SUBSET_ATTR]
        return cached
    setattr(context, _SUBSET_ATTR, None)
    vcf_samples = list(getattr(context, "vcf_samples", None) or [])
    if cfg.get("no_sample_subset") or not (cfg.get("perform_association") and vcf_samples):
        return None
    analysed = _analysed_samples(context, vcf_samples)
    if analysed is None:
        return None
    keep = [i for i, s in enumerate(vcf_samples) if s in analysed]
    if len(keep) == len(vcf_samples):
        return None
    reason = blocking_reason(cfg)
    if reason:
        logger.info(
            "Sample subsetting at extraction skipped (%s): keeping all %d VCF samples, %d analysed",
            reason,
            len(vcf_samples),
            len(keep),
        )
        return None
    raw = vcf_samples
    if (cfg.get("remove_sample_substring") or "").strip():
        from .helpers import get_vcf_samples

        raw = get_vcf_samples(cfg["vcf_file"])  # bcftools needs the header names
    subset = [vcf_samples[i] for i in keep]
    setattr(context, _VCF_COUNT_ATTR, len(vcf_samples))
    context.vcf_samples = subset
    cfg["sample_list"] = ",".join(subset)
    cfg["sample_subset_sha256"] = sample_set_sha256(subset)
    raw_subset = [raw[i] for i in keep]
    setattr(context, _SUBSET_ATTR, raw_subset)
    logger.info(
        "Sample subsetting at extraction: %d of %d VCF samples (cases/controls or phenotyped; "
        "sha256 %s)",
        len(subset),
        len(vcf_samples),
        cfg["sample_subset_sha256"][:12],
    )
    return raw_subset


def vcf_sample_count(context: Any) -> int | None:
    """VCF sample count before subsetting, or None when no subset was applied."""
    count: int | None = vars(context).get(_VCF_COUNT_ATTR)
    return count


def subset_summary(context: Any) -> str:
    """One-line statement of the sample-subset decision, from what the run already recorded."""
    cfg = context.config
    if sha := cfg.get("sample_subset_sha256"):
        return (
            f"sample subset: {len(context.vcf_samples or [])} VCF samples extracted (sha {sha[:8]})"
        )
    if cfg.get("no_sample_subset"):
        why = "--no-sample-subset"
    else:
        why = blocking_reason(cfg) or "no sample file or every VCF sample analysed"
    return f"sample subset: not used - {why}"
