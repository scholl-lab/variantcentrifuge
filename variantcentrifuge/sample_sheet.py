"""Validated sample-sheet adapter consolidating pedigree, phenotype, and BAM metadata."""

from __future__ import annotations

import argparse
import csv
import logging
import os
import re
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

logger = logging.getLogger("variantcentrifuge")

HPO_TERM_RE = re.compile(r"^HP:\d{7}$")

_SEX_MAP: dict[str, str] = {
    "": "0",
    ".": "0",
    "0": "0",
    "u": "0",
    "unk": "0",
    "unknown": "0",
    "1": "1",
    "m": "1",
    "male": "1",
    "2": "2",
    "f": "2",
    "female": "2",
}

_AFFECTED_MAP: dict[str, str] = {
    "": "0",
    ".": "0",
    "0": "0",
    "-9": "0",
    "u": "0",
    "unk": "0",
    "unknown": "0",
    "1": "1",
    "unaffected": "1",
    "control": "1",
    "2": "2",
    "affected": "2",
    "case": "2",
}


@dataclass
class SampleSheetData:
    """Consolidated sample metadata parsed from a validated sample sheet."""

    sample_ids: list[str] = field(default_factory=list)
    pedigree_data: dict[str, dict[str, Any]] = field(default_factory=dict)
    phenotype_data: dict[str, set[str]] = field(default_factory=dict)
    bam_mapping: dict[str, str] = field(default_factory=dict)


def _normalize_sex(raw_sex: str, sample_id: str) -> str:
    val = (raw_sex or "").strip()
    mapped = _SEX_MAP.get(val.lower())
    if mapped is None:
        raise ValueError(
            f"Invalid sex value '{val}' for sample '{sample_id}' in sample sheet "
            "(expected 0/1/2 or male/female/unknown)."
        )
    return mapped


def _normalize_affected_status(raw_status: str, sample_id: str) -> str:
    val = (raw_status or "").strip()
    mapped = _AFFECTED_MAP.get(val.lower())
    if mapped is None:
        raise ValueError(
            f"Invalid affected_status value '{val}' for sample '{sample_id}' in sample sheet "
            "(expected 0/1/2 or unaffected/affected/unknown)."
        )
    return mapped


def _parse_hpo_terms(raw_terms: str, sample_id: str) -> set[str]:
    terms: set[str] = set()
    if not raw_terms or not raw_terms.strip():
        return terms
    for token in re.split(r"[;,]", raw_terms):
        cleaned = token.strip()
        if not cleaned or cleaned.lower() in (".", "na", "none", "null"):
            continue
        if not HPO_TERM_RE.match(cleaned):
            raise ValueError(
                f"Malformed HPO term '{cleaned}' for sample '{sample_id}' in sample sheet "
                "(expected HP:NNNNNNN)."
            )
        terms.add(cleaned)
    return terms


def load_sample_sheet(file_path: str | Path) -> SampleSheetData:
    """Load and validate a TSV or CSV sample sheet."""
    path = Path(file_path)
    if not path.is_file():
        raise FileNotFoundError(f"Sample sheet file not found: {path}")
    if path.stat().st_size == 0:
        raise ValueError(f"Sample sheet file is empty: {path}")

    with open(path, encoding="utf-8") as handle:
        non_comment_lines = [
            line for line in handle if line.strip() and not line.lstrip().startswith("#")
        ]

    if not non_comment_lines:
        raise ValueError(f"Sample sheet file contains no header or data: {path}")

    header_line = non_comment_lines[0]
    delimiter = "\t" if "\t" in header_line else ("," if "," in header_line else "\t")
    reader = csv.DictReader(non_comment_lines, delimiter=delimiter)
    fieldnames = [f.strip() for f in (reader.fieldnames or []) if f is not None]

    if "sample_id" not in fieldnames:
        raise ValueError(f"Sample sheet must contain a 'sample_id' column: {path}")

    sample_ids: list[str] = []
    seen_samples: set[str] = set()
    pedigree_data: dict[str, dict[str, Any]] = {}
    raw_phenotypes: dict[str, set[str]] = {}
    bam_mapping: dict[str, str] = {}
    any_hpo_present = False

    for row_idx, raw_row in enumerate(reader, start=2):
        row = {(k.strip() if k else ""): (v.strip() if v else "") for k, v in raw_row.items()}
        sample_id = row.get("sample_id", "").strip()
        if not sample_id:
            raise ValueError(f"Empty sample_id at line {row_idx} in sample sheet: {path}")
        if sample_id in seen_samples:
            raise ValueError(f"Duplicate sample_id '{sample_id}' in sample sheet: {path}")

        seen_samples.add(sample_id)
        sample_ids.append(sample_id)

        family_id = row.get("family_id", "").strip()
        if not family_id or family_id == ".":
            family_id = "FAM1"

        father_id = (row.get("paternal_id") or row.get("father_id") or "").strip()
        if not father_id or father_id == ".":
            father_id = "0"

        mother_id = (row.get("maternal_id") or row.get("mother_id") or "").strip()
        if not mother_id or mother_id == ".":
            mother_id = "0"

        sex = _normalize_sex(row.get("sex", ""), sample_id)
        affected_raw = row.get("affected_status") or row.get("phenotype_status") or ""
        affected_status = _normalize_affected_status(affected_raw, sample_id)

        pedigree_data[sample_id] = {
            "family_id": family_id,
            "sample_id": sample_id,
            "father_id": father_id,
            "mother_id": mother_id,
            "sex": sex,
            "affected_status": affected_status,
        }

        hpo_raw = row.get("hpo_terms") or row.get("phenotypes") or row.get("phenotype") or ""
        parsed_terms = _parse_hpo_terms(hpo_raw, sample_id)
        if parsed_terms:
            any_hpo_present = True
        raw_phenotypes[sample_id] = parsed_terms

        bam_val = (row.get("bam_path") or row.get("bam") or "").strip()
        if bam_val and bam_val != ".":
            bam_mapping[sample_id] = bam_val

    if not sample_ids:
        raise ValueError(f"Sample sheet contains header only and no sample rows: {path}")

    phenotype_data = raw_phenotypes if any_hpo_present else {}
    return SampleSheetData(
        sample_ids=sample_ids,
        pedigree_data=pedigree_data,
        phenotype_data=phenotype_data,
        bam_mapping=bam_mapping,
    )


def validate_sample_sheet_conflicts(
    args: argparse.Namespace,
    cfg: dict[str, Any],
    parser: argparse.ArgumentParser | None = None,
) -> None:
    """Validate --sample-sheet and ensure it does not conflict with legacy metadata flags."""
    sample_sheet = getattr(args, "sample_sheet", None) or cfg.get("sample_sheet")
    if not sample_sheet:
        return

    has_ped = bool(
        getattr(args, "ped", None)
        or getattr(args, "ped_file", None)
        or cfg.get("ped_file")
        or cfg.get("ped")
    )
    has_pheno = bool(getattr(args, "phenotype_file", None) or cfg.get("phenotype_file"))
    has_bam = bool(
        (getattr(args, "bam_mapping_file", None) or cfg.get("bam_mapping_file"))
        and not cfg.get("_bam_mapping_from_sample_sheet")
    )

    if has_ped or has_pheno or has_bam:
        msg = (
            "--sample-sheet cannot be combined with --ped, --phenotype-file, or --bam-mapping-file"
        )
        if parser is not None:
            raise ValueError(msg)
        raise ValueError(msg)

    load_sample_sheet(str(sample_sheet))


def materialize_sample_sheet_bam_mapping(
    sheet: SampleSheetData,
    output_dir: str | Path,
) -> str | None:
    """Write a 2-column BAM mapping TSV in output_dir if the sample sheet contains BAM paths."""
    if not sheet.bam_mapping:
        return None
    out_dir = Path(output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    target = out_dir / "sample_sheet_bam_mapping.tsv"

    fd, tmp_name = tempfile.mkstemp(
        prefix=".sample_sheet_bam_mapping.",
        suffix=".tmp",
        dir=str(out_dir),
    )
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as handle:
            for sample_id, bam_path in sheet.bam_mapping.items():
                handle.write(f"{sample_id}\t{bam_path}\n")
        os.replace(tmp_name, target)
    except Exception:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)
        raise

    return str(target)


def classify_samples_from_phenotype_data(
    vcf_samples: list[str],
    phenotype_data: dict[str, set[str]],
    case_phenotypes: list[str],
    control_phenotypes: list[str],
    remove_substring: str = "",
) -> tuple[list[str], list[str]]:
    """Classify VCF samples into case and control lists from loaded phenotype_data."""
    case_phenotype_set = set(case_phenotypes) if case_phenotypes else set()
    control_phenotype_set = set(control_phenotypes) if control_phenotypes else set()

    vcf_samples_str = {str(s) for s in vcf_samples}
    if remove_substring and remove_substring.strip():
        vcf_samples_str = {s.replace(remove_substring, "") for s in vcf_samples_str}

    case_samples: set[str] = set()
    control_samples: set[str] = set()

    for raw_sample, sample_terms in phenotype_data.items():
        sample_key = str(raw_sample)
        if remove_substring and remove_substring.strip():
            sample_key = sample_key.replace(remove_substring, "")
        if case_phenotype_set and (sample_terms & case_phenotype_set):
            case_samples.add(sample_key)
        if control_phenotype_set and (sample_terms & control_phenotype_set):
            control_samples.add(sample_key)

    if case_phenotype_set and not control_phenotype_set:
        final_case_samples = case_samples & vcf_samples_str
        final_control_samples = vcf_samples_str - final_case_samples
    elif control_phenotype_set and not case_phenotype_set:
        final_control_samples = control_samples & vcf_samples_str
        final_case_samples = vcf_samples_str - final_control_samples
    else:
        final_case_samples = case_samples & vcf_samples_str
        final_control_samples = control_samples & vcf_samples_str

    return list(final_case_samples), list(final_control_samples)
