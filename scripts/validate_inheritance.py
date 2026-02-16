#!/usr/bin/env python3
"""
Golden file validation for inheritance analysis vectorization.

This script generates golden reference files from the current (pre-vectorization)
implementation and compares new outputs against them to ensure clinical equivalence.

Usage:
    python scripts/validate_inheritance.py generate
    python scripts/validate_inheritance.py compare
"""

import argparse
import json
import logging
import sys
from pathlib import Path

import pandas as pd

# Add project root to path
sys.path.insert(0, str(Path(__file__).parent.parent))

from variantcentrifuge.inheritance.analyzer import analyze_inheritance

logging.basicConfig(level=logging.INFO, format="%(levelname)s: %(message)s")
logger = logging.getLogger(__name__)

GOLDEN_DIR = Path(__file__).parent.parent / "tests" / "fixtures" / "golden"


def build_trio_denovo_scenario() -> tuple[pd.DataFrame, dict, list, str]:
    """Build trio with de novo variant (child het, parents ref)."""
    df = pd.DataFrame(
        [
            {
                "CHROM": "1",
                "POS": 1000,
                "REF": "A",
                "ALT": "T",
                "GENE": "GENE1",
                "child": "0/1",
                "father": "0/0",
                "mother": "0/0",
            }
        ]
    )

    pedigree = {
        "child": {
            "sample_id": "child",
            "father_id": "father",
            "mother_id": "mother",
            "sex": "1",
            "affected_status": "2",
        },
        "father": {
            "sample_id": "father",
            "father_id": "0",
            "mother_id": "0",
            "sex": "1",
            "affected_status": "1",
        },
        "mother": {
            "sample_id": "mother",
            "father_id": "0",
            "mother_id": "0",
            "sex": "2",
            "affected_status": "1",
        },
    }

    samples = ["child", "father", "mother"]
    return df, pedigree, samples, "trio_denovo"


def build_trio_dominant_scenario() -> tuple[pd.DataFrame, dict, list, str]:
    """Build trio with autosomal dominant (child het, father het+affected)."""
    df = pd.DataFrame(
        [
            {
                "CHROM": "2",
                "POS": 2000,
                "REF": "C",
                "ALT": "G",
                "GENE": "GENE2",
                "child": "0/1",
                "father": "0/1",
                "mother": "0/0",
            }
        ]
    )

    pedigree = {
        "child": {
            "sample_id": "child",
            "father_id": "father",
            "mother_id": "mother",
            "sex": "1",
            "affected_status": "2",
        },
        "father": {
            "sample_id": "father",
            "father_id": "0",
            "mother_id": "0",
            "sex": "1",
            "affected_status": "2",  # Father affected
        },
        "mother": {
            "sample_id": "mother",
            "father_id": "0",
            "mother_id": "0",
            "sex": "2",
            "affected_status": "1",
        },
    }

    samples = ["child", "father", "mother"]
    return df, pedigree, samples, "trio_dominant"


def build_trio_recessive_scenario() -> tuple[pd.DataFrame, dict, list, str]:
    """Build trio with autosomal recessive (child hom_alt, both parents het)."""
    df = pd.DataFrame(
        [
            {
                "CHROM": "3",
                "POS": 3000,
                "REF": "G",
                "ALT": "A",
                "GENE": "GENE3",
                "child": "1/1",
                "father": "0/1",
                "mother": "0/1",
            }
        ]
    )

    pedigree = {
        "child": {
            "sample_id": "child",
            "father_id": "father",
            "mother_id": "mother",
            "sex": "1",
            "affected_status": "2",
        },
        "father": {
            "sample_id": "father",
            "father_id": "0",
            "mother_id": "0",
            "sex": "1",
            "affected_status": "1",
        },
        "mother": {
            "sample_id": "mother",
            "father_id": "0",
            "mother_id": "0",
            "sex": "2",
            "affected_status": "1",
        },
    }

    samples = ["child", "father", "mother"]
    return df, pedigree, samples, "trio_recessive"


def build_trio_denovo_candidate_scenario() -> tuple[pd.DataFrame, dict, list, str]:
    """Build trio with de novo candidate (one parent missing GT)."""
    df = pd.DataFrame(
        [
            {
                "CHROM": "4",
                "POS": 4000,
                "REF": "T",
                "ALT": "C",
                "GENE": "GENE4",
                "child": "0/1",
                "father": "./.",
                "mother": "0/0",
            }
        ]
    )

    pedigree = {
        "child": {
            "sample_id": "child",
            "father_id": "father",
            "mother_id": "mother",
            "sex": "1",
            "affected_status": "2",
        },
        "father": {
            "sample_id": "father",
            "father_id": "0",
            "mother_id": "0",
            "sex": "1",
            "affected_status": "1",
        },
        "mother": {
            "sample_id": "mother",
            "father_id": "0",
            "mother_id": "0",
            "sex": "2",
            "affected_status": "1",
        },
    }

    samples = ["child", "father", "mother"]
    return df, pedigree, samples, "trio_denovo_candidate"


def build_single_sample_scenario() -> tuple[pd.DataFrame, dict, list, str]:
    """Build single sample (no pedigree): het -> unknown, hom_alt -> homozygous."""
    df = pd.DataFrame(
        [
            {
                "CHROM": "5",
                "POS": 5000,
                "REF": "A",
                "ALT": "G",
                "GENE": "GENE5",
                "Sample1": "0/1",
            },
            {
                "CHROM": "5",
                "POS": 5500,
                "REF": "C",
                "ALT": "T",
                "GENE": "GENE5",
                "Sample1": "1/1",
            },
        ]
    )

    pedigree = {
        "Sample1": {
            "sample_id": "Sample1",
            "father_id": "0",
            "mother_id": "0",
            "sex": "0",
            "affected_status": "2",
        }
    }

    samples = ["Sample1"]
    return df, pedigree, samples, "single_sample"


def build_extended_family_scenario() -> tuple[pd.DataFrame, dict, list, str]:
    """Build extended family (4+ members): dominant segregation across 3 generations."""
    df = pd.DataFrame(
        [
            {
                "CHROM": "6",
                "POS": 6000,
                "REF": "G",
                "ALT": "T",
                "GENE": "GENE6",
                "grandparent": "0/1",
                "parent": "0/1",
                "grandparent_spouse": "0/0",
                "child": "0/1",
                "spouse": "0/0",
            }
        ]
    )

    pedigree = {
        "grandparent": {
            "sample_id": "grandparent",
            "father_id": "0",
            "mother_id": "0",
            "sex": "1",
            "affected_status": "2",
        },
        "parent": {
            "sample_id": "parent",
            "father_id": "grandparent",
            "mother_id": "grandparent_spouse",
            "sex": "2",
            "affected_status": "2",
        },
        "grandparent_spouse": {
            "sample_id": "grandparent_spouse",
            "father_id": "0",
            "mother_id": "0",
            "sex": "2",
            "affected_status": "1",
        },
        "child": {
            "sample_id": "child",
            "father_id": "spouse",
            "mother_id": "parent",
            "sex": "1",
            "affected_status": "2",
        },
        "spouse": {
            "sample_id": "spouse",
            "father_id": "0",
            "mother_id": "0",
            "sex": "1",
            "affected_status": "1",
        },
    }

    samples = ["grandparent", "parent", "grandparent_spouse", "child", "spouse"]
    return df, pedigree, samples, "extended_family"


def build_x_linked_scenario() -> tuple[pd.DataFrame, dict, list, str]:
    """Build X-linked variants: CHROM=X, male proband het (hemizygous), mother carrier."""
    df = pd.DataFrame(
        [
            {
                "CHROM": "X",
                "POS": 7000,
                "REF": "A",
                "ALT": "C",
                "GENE": "DMD",
                "son": "0/1",
                "father": "0/0",
                "mother": "0/1",
            }
        ]
    )

    pedigree = {
        "son": {
            "sample_id": "son",
            "father_id": "father",
            "mother_id": "mother",
            "sex": "1",  # Male
            "affected_status": "2",
        },
        "father": {
            "sample_id": "father",
            "father_id": "0",
            "mother_id": "0",
            "sex": "1",
            "affected_status": "1",
        },
        "mother": {
            "sample_id": "mother",
            "father_id": "0",
            "mother_id": "0",
            "sex": "2",
            "affected_status": "1",  # Carrier
        },
    }

    samples = ["son", "father", "mother"]
    return df, pedigree, samples, "x_linked"


def build_mitochondrial_scenario() -> tuple[pd.DataFrame, dict, list, str]:
    """Build mitochondrial variants: CHROM=MT, maternal transmission."""
    df = pd.DataFrame(
        [
            {
                "CHROM": "MT",
                "POS": 8000,
                "REF": "T",
                "ALT": "G",
                "GENE": "MT-ND1",
                "child": "1/1",
                "father": "0/0",
                "mother": "1/1",
            }
        ]
    )

    pedigree = {
        "child": {
            "sample_id": "child",
            "father_id": "father",
            "mother_id": "mother",
            "sex": "1",
            "affected_status": "2",
        },
        "father": {
            "sample_id": "father",
            "father_id": "0",
            "mother_id": "0",
            "sex": "1",
            "affected_status": "1",
        },
        "mother": {
            "sample_id": "mother",
            "father_id": "0",
            "mother_id": "0",
            "sex": "2",
            "affected_status": "2",
        },
    }

    samples = ["child", "father", "mother"]
    return df, pedigree, samples, "mitochondrial"


def build_compound_het_scenario() -> tuple[pd.DataFrame, dict, list, str]:
    """Build compound het: Two het variants in same gene, one from each parent (trans config)."""
    df = pd.DataFrame(
        [
            {
                "CHROM": "7",
                "POS": 9000,
                "REF": "A",
                "ALT": "T",
                "GENE": "GENE7",
                "child": "0/1",
                "father": "0/1",
                "mother": "0/0",
            },
            {
                "CHROM": "7",
                "POS": 9500,
                "REF": "C",
                "ALT": "G",
                "GENE": "GENE7",
                "child": "0/1",
                "father": "0/0",
                "mother": "0/1",
            },
        ]
    )

    pedigree = {
        "child": {
            "sample_id": "child",
            "father_id": "father",
            "mother_id": "mother",
            "sex": "1",
            "affected_status": "2",
        },
        "father": {
            "sample_id": "father",
            "father_id": "0",
            "mother_id": "0",
            "sex": "1",
            "affected_status": "1",
        },
        "mother": {
            "sample_id": "mother",
            "father_id": "0",
            "mother_id": "0",
            "sex": "2",
            "affected_status": "1",
        },
    }

    samples = ["child", "father", "mother"]
    return df, pedigree, samples, "compound_het"


def build_edge_cases_scenario() -> tuple[pd.DataFrame, dict, list, str]:
    """Build edge cases: Missing genotypes, empty pedigree elements, single variant per gene."""
    df = pd.DataFrame(
        [
            {
                "CHROM": "8",
                "POS": 10000,
                "REF": "G",
                "ALT": "A",
                "GENE": "GENE8",
                "child": "0/1",
                "father": "./.",
                "mother": "./.",
            },
            {
                "CHROM": "9",
                "POS": 11000,
                "REF": "T",
                "ALT": "C",
                "GENE": "GENE9",
                "child": "0/1",
                "father": "0/0",
                "mother": "0/0",
            },
        ]
    )

    pedigree = {
        "child": {
            "sample_id": "child",
            "father_id": "father",
            "mother_id": "mother",
            "sex": "1",
            "affected_status": "2",
        },
        "father": {
            "sample_id": "father",
            "father_id": "0",
            "mother_id": "0",
            "sex": "1",
            "affected_status": "1",
        },
        "mother": {
            "sample_id": "mother",
            "father_id": "0",
            "mother_id": "0",
            "sex": "2",
            "affected_status": "1",
        },
    }

    samples = ["child", "father", "mother"]
    return df, pedigree, samples, "edge_cases"


def get_all_scenarios() -> list[tuple[pd.DataFrame, dict, list, str]]:
    """Get all test scenarios."""
    return [
        build_trio_denovo_scenario(),
        build_trio_dominant_scenario(),
        build_trio_recessive_scenario(),
        build_trio_denovo_candidate_scenario(),
        build_single_sample_scenario(),
        build_extended_family_scenario(),
        build_x_linked_scenario(),
        build_mitochondrial_scenario(),
        build_compound_het_scenario(),
        build_edge_cases_scenario(),
    ]


def generate_golden_files() -> int:
    """Generate golden reference files from current implementation."""
    logger.info("Generating golden reference files...")

    # Create golden directory
    GOLDEN_DIR.mkdir(parents=True, exist_ok=True)

    scenarios = get_all_scenarios()
    success_count = 0

    for df, pedigree_data, sample_list, scenario_name in scenarios:
        try:
            logger.info(f"Processing scenario: {scenario_name}")

            # Run analyze_inheritance with current code
            result_df = analyze_inheritance(df, pedigree_data, sample_list)

            # Save as parquet (preserves dtypes exactly)
            parquet_path = GOLDEN_DIR / f"{scenario_name}.parquet"
            result_df.to_parquet(parquet_path, index=False)
            logger.info(f"  Saved parquet: {parquet_path}")

            # Also save human-readable TSV for inspection
            tsv_path = GOLDEN_DIR / f"{scenario_name}.tsv"
            result_df.to_csv(tsv_path, sep="\t", index=False)
            logger.info(f"  Saved TSV: {tsv_path}")

            success_count += 1

        except Exception as e:
            logger.error(f"  Failed to process {scenario_name}: {e}")
            import traceback

            traceback.print_exc()

    logger.info(f"\nGenerated {success_count}/{len(scenarios)} golden files")
    return 0 if success_count == len(scenarios) else 1


def compare_dataframes(
    golden_df: pd.DataFrame, new_df: pd.DataFrame, scenario_name: str
) -> tuple[bool, list[str]]:
    """Compare two DataFrames for inheritance equivalence."""
    errors = []

    # Sort both by stable key (CHROM, POS, REF, ALT)
    sort_cols = ["CHROM", "POS", "REF", "ALT"]
    golden_sorted = golden_df.sort_values(sort_cols).reset_index(drop=True)
    new_sorted = new_df.sort_values(sort_cols).reset_index(drop=True)

    # Check row count
    if len(golden_sorted) != len(new_sorted):
        errors.append(
            f"Row count mismatch: golden={len(golden_sorted)}, new={len(new_sorted)}"
        )
        return False, errors

    # Compare Inheritance_Pattern column (must match exactly)
    if "Inheritance_Pattern" not in golden_sorted.columns:
        errors.append("Golden file missing Inheritance_Pattern column")
        return False, errors

    if "Inheritance_Pattern" not in new_sorted.columns:
        errors.append("New output missing Inheritance_Pattern column")
        return False, errors

    pattern_mismatches = 0
    for idx in range(len(golden_sorted)):
        golden_pattern = golden_sorted.loc[idx, "Inheritance_Pattern"]
        new_pattern = new_sorted.loc[idx, "Inheritance_Pattern"]

        if golden_pattern != new_pattern:
            pattern_mismatches += 1
            variant_key = (
                f"{golden_sorted.loc[idx, 'CHROM']}:{golden_sorted.loc[idx, 'POS']}"
                f":{golden_sorted.loc[idx, 'REF']}>{golden_sorted.loc[idx, 'ALT']}"
            )
            errors.append(
                f"  Row {idx} ({variant_key}): Pattern mismatch - "
                f"golden={golden_pattern}, new={new_pattern}"
            )

    # Compare Inheritance_Details column (parse JSON and compare)
    if "Inheritance_Details" not in golden_sorted.columns:
        errors.append("Golden file missing Inheritance_Details column")
        return False, errors

    if "Inheritance_Details" not in new_sorted.columns:
        errors.append("New output missing Inheritance_Details column")
        return False, errors

    details_mismatches = 0
    for idx in range(len(golden_sorted)):
        try:
            golden_details = json.loads(golden_sorted.loc[idx, "Inheritance_Details"])
            new_details = json.loads(new_sorted.loc[idx, "Inheritance_Details"])

            variant_key = (
                f"{golden_sorted.loc[idx, 'CHROM']}:{golden_sorted.loc[idx, 'POS']}"
                f":{golden_sorted.loc[idx, 'REF']}>{golden_sorted.loc[idx, 'ALT']}"
            )

            # Compare primary_pattern
            if golden_details.get("primary_pattern") != new_details.get("primary_pattern"):
                details_mismatches += 1
                errors.append(
                    f"  Row {idx} ({variant_key}): primary_pattern mismatch - "
                    f"golden={golden_details.get('primary_pattern')}, "
                    f"new={new_details.get('primary_pattern')}"
                )

            # Compare all_patterns (sets should be equal)
            golden_patterns = set(golden_details.get("all_patterns", []))
            new_patterns = set(new_details.get("all_patterns", []))
            if golden_patterns != new_patterns:
                details_mismatches += 1
                errors.append(
                    f"  Row {idx} ({variant_key}): all_patterns mismatch - "
                    f"golden={golden_patterns}, new={new_patterns}"
                )

            # Compare confidence (within epsilon=0.001)
            golden_conf = golden_details.get("confidence", 0.0)
            new_conf = new_details.get("confidence", 0.0)
            if abs(golden_conf - new_conf) > 0.001:
                details_mismatches += 1
                errors.append(
                    f"  Row {idx} ({variant_key}): confidence mismatch - "
                    f"golden={golden_conf:.3f}, new={new_conf:.3f}"
                )

            # Compare samples_with_pattern (sample IDs should match)
            golden_samples = {
                s["sample_id"] for s in golden_details.get("samples_with_pattern", [])
            }
            new_samples = {
                s["sample_id"] for s in new_details.get("samples_with_pattern", [])
            }
            if golden_samples != new_samples:
                details_mismatches += 1
                errors.append(
                    f"  Row {idx} ({variant_key}): samples_with_pattern mismatch - "
                    f"golden={golden_samples}, new={new_samples}"
                )

        except json.JSONDecodeError as e:
            errors.append(f"  Row {idx}: JSON decode error - {e}")
            details_mismatches += 1

    if pattern_mismatches > 0 or details_mismatches > 0:
        errors.insert(
            0,
            f"Found {pattern_mismatches} pattern mismatches and "
            f"{details_mismatches} detail mismatches",
        )
        return False, errors

    return True, []


def compare_outputs() -> int:
    """Compare new output against golden files."""
    logger.info("Comparing outputs against golden files...")

    if not GOLDEN_DIR.exists():
        logger.error(f"Golden directory does not exist: {GOLDEN_DIR}")
        logger.error("Run 'python scripts/validate_inheritance.py generate' first")
        return 1

    scenarios = get_all_scenarios()
    all_passed = True
    results = []

    for df, pedigree_data, sample_list, scenario_name in scenarios:
        golden_path = GOLDEN_DIR / f"{scenario_name}.parquet"

        if not golden_path.exists():
            logger.error(f"Golden file missing: {golden_path}")
            all_passed = False
            results.append((scenario_name, False, [f"Golden file not found: {golden_path}"]))
            continue

        try:
            # Load golden file
            golden_df = pd.read_parquet(golden_path)

            # Run analyze_inheritance with current code
            new_df = analyze_inheritance(df, pedigree_data, sample_list)

            # Compare
            passed, errors = compare_dataframes(golden_df, new_df, scenario_name)
            results.append((scenario_name, passed, errors))

            if not passed:
                all_passed = False

        except Exception as e:
            logger.error(f"Error comparing {scenario_name}: {e}")
            import traceback

            traceback.print_exc()
            all_passed = False
            results.append((scenario_name, False, [f"Exception: {e}"]))

    # Print summary
    logger.info("\n" + "=" * 80)
    logger.info("COMPARISON RESULTS")
    logger.info("=" * 80)

    passed_count = sum(1 for _, passed, _ in results if passed)
    total_count = len(results)

    for scenario_name, passed, errors in results:
        status = "PASS" if passed else "FAIL"
        logger.info(f"{status:6s} - {scenario_name}")
        if errors:
            for error in errors:
                logger.info(f"         {error}")

    logger.info("=" * 80)
    logger.info(f"TOTAL: {passed_count}/{total_count} scenarios passed")

    return 0 if all_passed else 1


def main() -> int:
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Golden file validation for inheritance analysis"
    )
    parser.add_argument(
        "mode",
        choices=["generate", "compare"],
        help="Mode: generate golden files or compare against them",
    )

    args = parser.parse_args()

    if args.mode == "generate":
        return generate_golden_files()
    elif args.mode == "compare":
        return compare_outputs()
    else:
        parser.print_help()
        return 1


if __name__ == "__main__":
    sys.exit(main())
