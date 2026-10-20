# File: variantcentrifuge/association/meta/cli.py
"""
``variantcentrifuge meta`` subcommand.

Consumes pre-computed per-stratum score/covariance artifacts (no VCF, no external
tools) and writes a meta-analysis results TSV. This mirrors the RAREMETAL /
MetaSKAT / seqMeta pattern: per-study score files + a separate meta step.
"""

from __future__ import annotations

import argparse
import logging
import sys

logger = logging.getLogger("variantcentrifuge")


def build_parser() -> argparse.ArgumentParser:
    """Argument parser for the meta-analysis CLI."""
    parser = argparse.ArgumentParser(
        prog="variantcentrifuge meta",
        description=(
            "Cohort-stratified rare-variant score-statistic meta-analysis. "
            "Combines per-stratum score/covariance artifacts (binary traits)."
        ),
    )
    parser.add_argument(
        "--manifest", required=True, help="Study manifest TSV (stratum_id, artifact_path)."
    )
    parser.add_argument("--output", required=True, help="Output meta results TSV path.")
    parser.add_argument(
        "--weights",
        default="beta:1,25",
        help="Weight scheme: 'beta:a,b' (default beta:1,25) or 'uniform'.",
    )
    parser.add_argument(
        "--correction",
        default="fdr",
        choices=["fdr", "bonferroni"],
        help="Multiple-testing correction across genes (default: fdr).",
    )
    parser.add_argument(
        "--min-mac",
        type=float,
        default=0.0,
        help="Genes with cumulative minor-allele count below N get status=low_mac and "
        "are excluded from the multiple-testing correction (default: 0, no filter).",
    )
    parser.add_argument(
        "--require-variants-in-all-strata",
        action="store_true",
        help="Genes without an informative minor allele in every stratum get "
        "status=not_in_all_strata and are excluded from the correction.",
    )
    parser.add_argument(
        "--diagnostics-output",
        default=None,
        metavar="DIR",
        help="Write lambda_GC and QQ diagnostics for the meta p-values to DIR.",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging verbosity (default: INFO).",
    )
    return parser


def main(argv: list[str] | None = None) -> int:
    """CLI entry point. Returns a process exit code."""
    args = build_parser().parse_args(argv)

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s [%(levelname)s] %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stderr,
    )

    # Import here so --help stays fast and import errors surface at run time.
    from variantcentrifuge.association.meta.manifest import MetaCompatibilityError
    from variantcentrifuge.association.meta.runner import run_meta

    try:
        df = run_meta(
            manifest_path=args.manifest,
            output_path=args.output,
            weight_scheme=args.weights,
            correction=args.correction,
            min_mac=args.min_mac,
            diagnostics_dir=args.diagnostics_output,
            require_variants_in_all_strata=args.require_variants_in_all_strata,
        )
    except MetaCompatibilityError as exc:
        logger.error("Meta-analysis refused (strata not comparable): %s", exc)
        return 2
    except (FileNotFoundError, ValueError) as exc:
        logger.error("Meta-analysis failed: %s", exc)
        return 1

    logger.info("Meta-analysis wrote %d gene rows to %s", len(df), args.output)
    return 0


if __name__ == "__main__":  # pragma: no cover
    raise SystemExit(main())
