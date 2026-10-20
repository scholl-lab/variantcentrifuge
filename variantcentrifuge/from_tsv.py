"""Support for running VariantCentrifuge downstream stages from a pre-extracted/final TSV."""

from __future__ import annotations

import argparse
import gzip
import logging
import os
from pathlib import Path
from typing import TYPE_CHECKING, Any, cast

import pandas as pd

from .gt_reconstruction import _find_per_sample_gt_columns
from .links import resolve_links
from .ped_reader import read_pedigree
from .phenotype import load_phenotypes
from .phenotype_helpers import extract_sample_and_genotype

if TYPE_CHECKING:
    from .pipeline_core import PipelineContext

logger = logging.getLogger("variantcentrifuge")

REQUIRED_BASE_COLUMNS = ("CHROM", "POS", "REF", "ALT")
ACCEPTED_GENE_COLUMNS = ("GENE", "ANN[0].GENE", "ANN_0__GENE")

EXTRACTED_TO_STANDARD_COLUMNS = {
    "ANN[0].GENE": "GENE",
    "ANN_0__GENE": "GENE",
    "ANN[0].EFFECT": "EFFECT",
    "ANN_0__EFFECT": "EFFECT",
    "ANN[0].IMPACT": "IMPACT",
    "ANN_0__IMPACT": "IMPACT",
    "ANN[0].HGVS_C": "HGVS_C",
    "ANN_0__HGVS_C": "HGVS_C",
    "ANN[0].HGVS_P": "HGVS_P",
    "ANN_0__HGVS_P": "HGVS_P",
    "ANN[0].FEATUREID": "FEATUREID",
    "ANN_0__FEATUREID": "FEATUREID",
}


def compute_tsv_base_name(tsv_file: str) -> str:
    """Derive a clean base name from a TSV or TSV.GZ path."""
    name = os.path.basename(tsv_file)
    for ext in (".tsv.gz", ".txt.gz", ".tsv", ".txt", ".gz"):
        if name.lower().endswith(ext):
            name = name[: -len(ext)]
            break
    return name or "from_tsv"


def _collect_metadata_samples(cfg: dict[str, Any]) -> list[str]:
    """Collect ordered unique sample IDs from PED, phenotype, or case/control config."""
    ordered: list[str] = []
    seen: set[str] = set()

    def _add(sample_id: str) -> None:
        s = sample_id.strip()
        if s and s not in seen:
            seen.add(s)
            ordered.append(s)

    sample_sheet = cfg.get("sample_sheet")
    if sample_sheet and os.path.isfile(str(sample_sheet)):
        from .sample_sheet import load_sample_sheet

        sheet_data = load_sample_sheet(str(sample_sheet))
        for sid in sheet_data.sample_ids:
            _add(str(sid))

    ped_file = cfg.get("ped_file") or cfg.get("ped")
    if ped_file and os.path.isfile(str(ped_file)):
        ped_data = read_pedigree(str(ped_file))
        for sid in ped_data:
            _add(str(sid))

    pheno_file = cfg.get("phenotype_file")
    pheno_sample_col = cfg.get("phenotype_sample_column")
    pheno_value_col = cfg.get("phenotype_value_column")
    if pheno_file and os.path.isfile(str(pheno_file)):
        if pheno_sample_col and pheno_value_col:
            pheno_map = load_phenotypes(
                str(pheno_file), str(pheno_sample_col), str(pheno_value_col)
            )
            for sid in pheno_map:
                _add(str(sid))
        else:
            pheno_df = pd.read_csv(str(pheno_file), sep=None, engine="python", dtype=str)
            if not pheno_df.empty:
                first_col = pheno_df.columns[0]
                for val in pheno_df[first_col].dropna():
                    _add(str(val))

    for list_key, file_key in (
        ("case_samples", "case_samples_file"),
        ("control_samples", "control_samples_file"),
        ("vcf_samples", None),
    ):
        samples_val = cfg.get(list_key)
        if isinstance(samples_val, str):
            samples_val = [s.strip() for s in samples_val.split(",") if s.strip()]
        if isinstance(samples_val, list):
            for sid in samples_val:
                _add(str(sid))
        file_val = cfg.get(file_key) if file_key else None
        if file_val and os.path.isfile(str(file_val)):
            with open(str(file_val), encoding="utf-8") as fh:
                for line in fh:
                    _add(line)

    return ordered


def _extract_gt_samples_from_series(gt_series: pd.Series) -> list[str]:
    """Extract ordered unique sample IDs from a packed GT column ('SampleA(0/1);SampleB(0/0)')."""
    ordered: list[str] = []
    seen: set[str] = set()
    for val in gt_series.dropna():
        text = str(val).strip()
        if not text:
            continue
        for entry in text.split(";"):
            entry = entry.strip()
            if not entry:
                continue
            sample_name, _gt = extract_sample_and_genotype(entry)
            if sample_name and sample_name not in seen:
                seen.add(sample_name)
                ordered.append(sample_name)
    return ordered


def infer_and_validate_from_tsv_samples(df: pd.DataFrame, cfg: dict[str, Any]) -> list[str]:
    """Infer and validate the ordered sample list from a TSV DataFrame and sample metadata."""
    meta_samples = _collect_metadata_samples(cfg)
    anon_gt_cols = _find_per_sample_gt_columns(df)

    if anon_gt_cols:
        n_gen = len(anon_gt_cols)
        if meta_samples:
            if len(meta_samples) != n_gen:
                raise ValueError(
                    f"Ambiguous sample mapping in --from-tsv mode: TSV contains {n_gen} "
                    f"anonymous per-sample GT columns ({anon_gt_cols[0]}..{anon_gt_cols[-1]}), "
                    f"but supplied sample metadata defines {len(meta_samples)} sample(s) "
                    f"({meta_samples}). Counts must match exactly."
                )
            return meta_samples

        needs_named_samples = bool(
            cfg.get("calculate_inheritance")
            or cfg.get("ped_file")
            or cfg.get("phenotype_file")
            or cfg.get("perform_gene_burden")
            or cfg.get("perform_association")
        )
        if needs_named_samples:
            raise ValueError(
                "Ambiguous sample mapping in --from-tsv mode: TSV contains anonymous "
                f"per-sample GT columns ({anon_gt_cols[0]}..), which require sample metadata "
                "(--ped, --phenotype-file, --sample-sheet, or --case-samples/--control-samples) "
                "to map columns to sample IDs."
            )
        return [f"SAMPLE_{i}" for i in range(n_gen)]

    if "GT" in df.columns:
        gt_samples = _extract_gt_samples_from_series(df["GT"])
        if gt_samples and meta_samples and not (set(gt_samples) & set(meta_samples)):
            raise ValueError(
                "Sample IDs in --from-tsv GT column do not overlap with supplied sample "
                f"metadata: GT samples={gt_samples[:5]}, metadata samples={meta_samples[:5]}."
            )
        if meta_samples:
            meta_set = set(meta_samples)
            return meta_samples + [s for s in gt_samples if s not in meta_set]
        if gt_samples:
            return gt_samples
        return ["SAMPLE_0"] if len(df) > 0 else []

    if meta_samples and any(s in df.columns for s in meta_samples):
        return [s for s in meta_samples if s in df.columns]

    raise ValueError(
        "Invalid --from-tsv file: missing genotype column ('GT' or per-sample 'GEN[i].GT')."
    )


def validate_from_tsv_args_and_file(
    args: argparse.Namespace,
    cfg: dict[str, Any],
    parser: argparse.ArgumentParser | None = None,
) -> list[str]:
    """Validate CLI arguments and input TSV file for --from-tsv mode."""
    from_tsv = getattr(args, "from_tsv", None) or cfg.get("from_tsv")
    if not from_tsv:
        raise ValueError("--from-tsv path is required in TSV mode.")

    if getattr(args, "vcf_file", None):
        msg = "Cannot specify both -v/--vcf-file and --from-tsv."
        if parser is not None:
            parser.error(msg)
        raise ValueError(msg)

    # Reject VCF-only flags
    gene_name = getattr(args, "gene_name", None)
    vcf_only_checks = [
        ("--bcftools-prefilter", getattr(args, "bcftools_prefilter", None)),
        ("--split-snpeff-lines", getattr(args, "split_snpeff_lines", None)),
        ("--transcript-list", getattr(args, "transcript_list", None)),
        ("--transcript-file", getattr(args, "transcript_file", None)),
        (
            "-g/--gene-name",
            None if not gene_name or str(gene_name).strip().lower() == "all" else gene_name,
        ),
        ("-G/--gene-file", getattr(args, "gene_file", None)),
        ("--regions-bed", getattr(args, "regions_bed", None)),
    ]
    for flag_name, val in vcf_only_checks:
        if val:
            raise ValueError(
                f"{flag_name} requires VCF extraction and cannot be used with --from-tsv."
            )

    if (getattr(args, "pca", None) or cfg.get("pca")) == "akt":
        raise ValueError("--pca akt requires a VCF file and cannot be used with --from-tsv.")

    tsv_path = Path(str(from_tsv))
    if not tsv_path.exists():
        raise FileNotFoundError(f"Input TSV file does not exist: {tsv_path}")
    if not tsv_path.is_file():
        raise ValueError(f"Input TSV path is not a regular file: {tsv_path}")
    if tsv_path.stat().st_size == 0:
        raise ValueError(f"Input TSV file is empty (0 bytes): {tsv_path}")

    open_fn = gzip.open if str(tsv_path).endswith(".gz") else open
    try:
        with open_fn(tsv_path, "rt", encoding="utf-8-sig") as fh:
            header_line = fh.readline().rstrip("\r\n")
    except Exception as e:
        raise ValueError(f"Failed to read --from-tsv file '{tsv_path}': {e}") from e

    if not header_line.strip():
        raise ValueError(f"Input TSV file has no header row: {tsv_path}")

    columns = header_line.split("\t")
    missing_base = [c for c in REQUIRED_BASE_COLUMNS if c not in columns]
    gene_col_cfg = cfg.get("gene_column_name")
    accepted_genes = (
        (*ACCEPTED_GENE_COLUMNS, str(gene_col_cfg)) if gene_col_cfg else ACCEPTED_GENE_COLUMNS
    )
    has_gene = any(c in columns for c in accepted_genes)
    if missing_base or not has_gene:
        missing_desc = list(missing_base)
        if not has_gene:
            missing_desc.append("GENE (or ANN[0].GENE)")
        raise ValueError(
            f"Input TSV '{tsv_path}' is missing required column(s): {', '.join(missing_desc)}"
        )

    compression = "gzip" if str(tsv_path).endswith(".gz") else None
    df = pd.read_csv(
        tsv_path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        encoding="utf-8-sig",
        compression=cast(Any, compression),
    )

    samples = infer_and_validate_from_tsv_samples(df, cfg)
    cfg["_from_tsv_samples"] = samples
    return samples


def prepare_from_tsv_context(context: PipelineContext) -> PipelineContext:
    """Prepare context.vcf_samples and normalized intermediate TSV for --from-tsv mode."""
    from_tsv = context.config.get("from_tsv") or getattr(context.args, "from_tsv", None)
    if not from_tsv:
        raise ValueError("No --from-tsv file configured in context.")

    tsv_path = Path(str(from_tsv))
    compression = "gzip" if str(tsv_path).endswith(".gz") else None
    df = pd.read_csv(
        tsv_path,
        sep="\t",
        dtype=str,
        keep_default_na=False,
        encoding="utf-8-sig",
        compression=cast(Any, compression),
    )

    # Normalize extracted column headers (ANN[0].GENE -> GENE, etc.) if short name not present
    rename_map: dict[str, str] = {}
    for src_col, dst_col in EXTRACTED_TO_STANDARD_COLUMNS.items():
        if (
            src_col in df.columns
            and dst_col not in df.columns
            and dst_col not in rename_map.values()
        ):
            rename_map[src_col] = dst_col
    if rename_map:
        df = df.rename(columns=rename_map)

    # Strip pre-existing external link columns so TSVOutputStage does not create .1 duplicates
    if not context.config.get("no_links", False):
        link_names = set(resolve_links(context.config).keys())
        cols_to_drop = [
            c
            for c in df.columns
            if c in link_names or any(c.startswith(f"{ln}.") for ln in link_names)
        ]
        if cols_to_drop:
            df = df.drop(columns=cols_to_drop)

    # Drop stale Phenotypes column when phenotype data is provided for recomputation
    if context.phenotype_data and "Phenotypes" in df.columns:
        df = df.drop(columns=["Phenotypes"])

    # Drop stale inheritance columns when inheritance calculation is requested
    if (
        context.config.get("calculate_inheritance")
        or context.config.get("ped_file")
        or context.config.get("sample_sheet")
        or context.config.get("inheritance_mode")
    ):
        inh_cols = [
            c
            for c in (
                "Inheritance_Pattern",
                "Inheritance_Details",
                "Inheritance_Confidence",
                "Inheritance_Description",
                "Inheritance_Samples",
            )
            if c in df.columns
        ]
        if inh_cols:
            df = df.drop(columns=inh_cols)

    samples = context.config.get("_from_tsv_samples")
    if samples is None:
        samples = infer_and_validate_from_tsv_samples(df, context.config)

    remove_substring = context.config.get("remove_sample_substring")
    if remove_substring and remove_substring.strip():
        samples = [s.replace(remove_substring, "") for s in samples]

    context.vcf_samples = list(samples)
    context.config["sample_list"] = ",".join(context.vcf_samples)

    prepared_tsv = context.workspace.get_intermediate_path(
        f"{context.workspace.base_name}.from_tsv_input.tsv"
    )
    if context.config.get("gzip_intermediates"):
        prepared_tsv = Path(str(prepared_tsv) + ".gz")

    out_compression = "gzip" if str(prepared_tsv).endswith(".gz") else None
    df.to_csv(prepared_tsv, sep="\t", index=False, compression=cast(Any, out_compression))

    context.extracted_tsv = prepared_tsv
    context.data = prepared_tsv
    return context
