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
    # VCF file and TSV input
    if getattr(args, "vcf_file", None):
        cfg["vcf_file"] = args.vcf_file
    elif cfg.get("vcf_file"):
        args.vcf_file = cfg["vcf_file"]

    if getattr(args, "from_tsv", None):
        cfg["from_tsv"] = args.from_tsv
    elif cfg.get("from_tsv"):
        args.from_tsv = cfg["from_tsv"]

    # Gene selection: explicit CLI flag clears opposing configured selector
    if getattr(args, "gene_name", None) is not None:
        cfg["gene_name"] = args.gene_name
        cfg["gene_file"] = None
        args.gene_file = None
    elif getattr(args, "gene_file", None) is not None:
        cfg["gene_file"] = args.gene_file
        cfg["gene_name"] = None
        args.gene_name = None
    else:
        if cfg.get("gene_name") is not None:
            args.gene_name = cfg["gene_name"]
        if cfg.get("gene_file") is not None:
            args.gene_file = cfg["gene_file"]

    # I/O, reporting, and execution flags
    for flag in (
        "xlsx",
        "html_report",
        "igv",
        "archive_results",
        "keep_intermediates",
        "late_filtering",
        "pseudonymize",
    ):
        flag_val = bool(getattr(args, flag, False) or cfg.get(flag, False))
        setattr(args, flag, flag_val)
        cfg[flag] = flag_val

    pseudonymize_table = getattr(args, "pseudonymize_table", None) or cfg.get("pseudonymize_table")
    if pseudonymize_table:
        args.pseudonymize_table = pseudonymize_table
        cfg["pseudonymize_table"] = pseudonymize_table
        args.pseudonymize = True
        cfg["pseudonymize"] = True

    gzip_inter = getattr(args, "gzip_intermediates", None)
    if gzip_inter is not None:
        cfg["gzip_intermediates"] = gzip_inter
    else:
        args.gzip_intermediates = cfg.get("gzip_intermediates", True)
        cfg["gzip_intermediates"] = args.gzip_intermediates

    sample_sheet = getattr(args, "sample_sheet", None) or cfg.get("sample_sheet")
    if sample_sheet:
        args.sample_sheet = sample_sheet
        cfg["sample_sheet"] = sample_sheet

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
        stat_flag_val = bool(getattr(args, flag, False) or cfg.get(flag, False))
        setattr(args, flag, stat_flag_val)
        cfg[flag] = stat_flag_val

    assoc_tests = getattr(args, "association_tests", None) or cfg.get("association_tests")
    if assoc_tests is not None:
        if isinstance(assoc_tests, str):
            assoc_tests = [t.strip() for t in assoc_tests.split(",") if t.strip()]
        args.association_tests = assoc_tests
        cfg["association_tests"] = assoc_tests

    # Inheritance analysis flags
    inh_mode = getattr(args, "inheritance_mode", None) or cfg.get("inheritance_mode")
    if inh_mode is not None:
        args.inheritance_mode = inh_mode
        cfg["inheritance_mode"] = inh_mode
    cfg["calculate_inheritance"] = bool(
        inh_mode is not None
        or ped is not None
        or sample_sheet is not None
        or cfg.get("calculate_inheritance", False)
    )
    args.calculate_inheritance = cfg["calculate_inheritance"]

    # PCA flags
    pca = getattr(args, "pca", None) or cfg.get("pca")
    if pca is not None:
        args.pca = pca
        cfg["pca"] = pca

    # Threads
    from .pipeline_validation import resolve_threads

    raw_threads = getattr(args, "threads", None) or cfg.get("threads", "auto")
    resolved_threads = resolve_threads(raw_threads, cfg)
    args.threads = resolved_threads
    cfg["threads"] = resolved_threads

    # Fields and IGV
    fields_file = getattr(args, "fields_file", None) or cfg.get("fields_file")
    if fields_file:
        args.fields_file = fields_file
        cfg["fields_file"] = fields_file

    # IGV companion parameters
    for opt in (
        "bam_mapping_file",
        "igv_reference",
        "igv_fasta",
        "igv_ideogram",
    ):
        igv_val = getattr(args, opt, None) or cfg.get(opt)
        if igv_val is not None:
            setattr(args, opt, igv_val)
            cfg[opt] = igv_val

    igv_flank = getattr(args, "igv_flanking", None)
    if igv_flank is not None:
        cfg["igv_flanking"] = int(igv_flank)
    elif "igv_flanking" in cfg:
        args.igv_flanking = int(cfg["igv_flanking"])

    max_depth = (
        getattr(args, "max_read_depth", None)
        or cfg.get("max_read_depth")
        or cfg.get("max_Read_depth")
    )
    if max_depth is not None:
        args.max_read_depth = max_depth
        cfg["max_read_depth"] = max_depth

    # Stage-activating filters and scoring options
    final_filter = getattr(args, "final_filter", None) or cfg.get("final_filter")
    if final_filter:
        args.final_filter = final_filter
        cfg["final_filter"] = final_filter

    tx_list = getattr(args, "transcript_list", None) or cfg.get("transcript_list")
    if tx_list:
        args.transcript_list = tx_list
        cfg["transcript_list"] = tx_list

    tx_file = getattr(args, "transcript_file", None) or cfg.get("transcript_file")
    if tx_file:
        args.transcript_file = tx_file
        cfg["transcript_file"] = tx_file

    clinvar = getattr(args, "clinvar_pm5_lookup", None) or cfg.get("clinvar_pm5_lookup")
    if clinvar:
        args.clinvar_pm5_lookup = clinvar
        cfg["clinvar_pm5_lookup"] = clinvar

    scoring_path = getattr(args, "scoring_config_path", None) or cfg.get("scoring_config_path")
    if scoring_path:
        args.scoring_config_path = scoring_path
        cfg["scoring_config_path"] = scoring_path

    # Custom annotation options
    ann_bed = (
        getattr(args, "annotate_bed", None)
        or cfg.get("annotate_bed")
        or cfg.get("annotate_bed_files")
    )
    if ann_bed:
        args.annotate_bed = ann_bed
        cfg["annotate_bed"] = ann_bed
        cfg["annotate_bed_files"] = ann_bed

    ann_gene = (
        getattr(args, "annotate_gene_list", None)
        or cfg.get("annotate_gene_list")
        or cfg.get("annotate_gene_lists")
    )
    if ann_gene:
        args.annotate_gene_list = ann_gene
        cfg["annotate_gene_list"] = ann_gene
        cfg["annotate_gene_lists"] = ann_gene

    ann_json = getattr(args, "annotate_json_genes", None) or cfg.get("annotate_json_genes")
    if ann_json:
        args.annotate_json_genes = ann_json
        cfg["annotate_json_genes"] = ann_json

    json_mapping = getattr(args, "json_gene_mapping", None) or cfg.get("json_gene_mapping")
    if json_mapping:
        args.json_gene_mapping = json_mapping
        cfg["json_gene_mapping"] = json_mapping

    # Stage discovery parity: multiallelic split, no-replacement, filters, presets
    split_snpeff = (
        getattr(args, "split_snpeff_lines", None)
        or cfg.get("split_snpeff_lines")
        or cfg.get("snpeff_splitting_mode")
    )
    if split_snpeff:
        args.split_snpeff_lines = split_snpeff
        cfg["split_snpeff_lines"] = split_snpeff
        cfg["snpeff_splitting_mode"] = split_snpeff
        cfg["snpeff_split_by_transcript"] = True

    no_rep = bool(getattr(args, "no_replacement", False) or cfg.get("no_replacement", False))
    args.no_replacement = no_rep
    cfg["no_replacement"] = no_rep

    filters = getattr(args, "filters", None) or cfg.get("filters")
    if filters:
        args.filters = filters
        cfg["filters"] = filters

    preset = getattr(args, "preset", None) or cfg.get("preset")
    if preset:
        args.preset = preset
        cfg["preset"] = preset

    for opt in (
        "remove_snpeff_sep",
        "do_not_remove_additional_columns",
        "extra_columns_to_remove",
    ):
        col_val = getattr(args, opt, None)
        if col_val is not None:
            cfg[opt] = col_val
        elif opt in cfg:
            setattr(args, opt, cfg[opt])
