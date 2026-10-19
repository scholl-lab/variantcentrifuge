"""CLI and configuration synchronization helpers."""

from __future__ import annotations

import argparse
from typing import Any


def sync_cli_and_config(args: argparse.Namespace, cfg: dict[str, Any]) -> None:
    """Synchronize settings between CLI arguments and configuration dictionary.

    Explicit CLI options override configuration file values, while configuration
    values populate omitted CLI options so downstream validation, execution, and
    stage discovery see consistent effective settings.
    """
    # VCF file
    if getattr(args, "vcf_file", None):
        cfg["vcf_file"] = args.vcf_file
    elif cfg.get("vcf_file"):
        args.vcf_file = cfg["vcf_file"]

    # Gene selection
    if getattr(args, "gene_name", None) is not None:
        cfg["gene_name"] = args.gene_name
    elif cfg.get("gene_name") is not None:
        args.gene_name = cfg["gene_name"]

    if getattr(args, "gene_file", None) is not None:
        cfg["gene_file"] = args.gene_file
    elif cfg.get("gene_file") is not None:
        args.gene_file = cfg["gene_file"]

    # I/O and reporting flags
    for flag in ("xlsx", "html_report", "archive_results", "keep_intermediates"):
        val = bool(getattr(args, flag, False) or cfg.get(flag, False))
        setattr(args, flag, val)
        cfg[flag] = val

    gzip_inter = getattr(args, "gzip_intermediates", None)
    if gzip_inter is not None:
        cfg["gzip_intermediates"] = gzip_inter
    else:
        args.gzip_intermediates = cfg.get("gzip_intermediates", True)
        cfg["gzip_intermediates"] = args.gzip_intermediates

    # Pedigree and phenotype files/columns
    ped = (
        getattr(args, "ped", None)
        or getattr(args, "ped_file", None)
        or cfg.get("ped_file")
        or cfg.get("ped")
    )
    if ped:
        args.ped = ped
        args.ped_file = ped
        cfg["ped_file"] = ped
        cfg["ped"] = ped

    pheno = getattr(args, "phenotype_file", None) or cfg.get("phenotype_file")
    if pheno:
        args.phenotype_file = pheno
        cfg["phenotype_file"] = pheno

    pheno_s_col = getattr(args, "phenotype_sample_column", None) or cfg.get(
        "phenotype_sample_column"
    )
    if pheno_s_col:
        args.phenotype_sample_column = pheno_s_col
        cfg["phenotype_sample_column"] = pheno_s_col

    pheno_v_col = getattr(args, "phenotype_value_column", None) or cfg.get("phenotype_value_column")
    if pheno_v_col:
        args.phenotype_value_column = pheno_v_col
        cfg["phenotype_value_column"] = pheno_v_col

    # Statistical / association analysis flags
    for flag in ("perform_gene_burden", "perform_association", "no_stats"):
        val = bool(getattr(args, flag, False) or cfg.get(flag, False))
        setattr(args, flag, val)
        cfg[flag] = val

    assoc_tests = getattr(args, "association_tests", None) or cfg.get("association_tests")
    if assoc_tests is not None:
        if isinstance(assoc_tests, str):
            assoc_tests = [t.strip() for t in assoc_tests.split(",") if t.strip()]
        args.association_tests = assoc_tests
        cfg["association_tests"] = assoc_tests

    # Threads
    threads = getattr(args, "threads", None) or cfg.get("threads", 1)
    args.threads = int(threads)
    cfg["threads"] = int(threads)

    # Fields and IGV
    fields_file = getattr(args, "fields_file", None) or cfg.get("fields_file")
    if fields_file:
        args.fields_file = fields_file
        cfg["fields_file"] = fields_file

    igv_flanking = getattr(args, "igv_flanking", None)
    if igv_flanking is not None:
        cfg["igv_flanking"] = int(igv_flanking)
    elif "igv_flanking" in cfg:
        args.igv_flanking = int(cfg["igv_flanking"])
