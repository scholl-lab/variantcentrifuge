#!/usr/bin/env python
"""Compare outputs between old and new pipeline implementations.

This script provides detailed comparison of pipeline outputs, helping
identify any differences between implementations.
"""

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd
from colorama import Fore, Style, init

# Initialize colorama for cross-platform colored output
init(autoreset=True)


class OutputComparator:
    """Compares outputs between pipeline implementations."""

    def __init__(self, verbose: bool = False):
        self.verbose = verbose

    def compare_directories(
        self,
        old_dir: Path,
        new_dir: Path,
        pattern: str = "*.tsv",
    ) -> Dict[str, Dict]:
        """Compare all matching files in two directories.

        Returns dict mapping filename to comparison results.
        """
        old_files = set(f.name for f in old_dir.glob(pattern))
        new_files = set(f.name for f in new_dir.glob(pattern))

        results = {}

        # Files only in old
        for f in old_files - new_files:
            results[f] = {"status": "missing_in_new", "details": []}

        # Files only in new
        for f in new_files - old_files:
            results[f] = {"status": "missing_in_old", "details": []}

        # Files in both
        for f in old_files & new_files:
            old_path = old_dir / f
            new_path = new_dir / f

            if f.endswith(".tsv"):
                match, details = self.compare_tsv_files(old_path, new_path)
                results[f] = {
                    "status": "match" if match else "differ",
                    "details": details,
                }
            elif f.endswith(".json"):
                match, details = self.compare_json_files(old_path, new_path)
                results[f] = {
                    "status": "match" if match else "differ",
                    "details": details,
                }
            else:
                # Binary comparison
                match = self._files_identical(old_path, new_path)
                results[f] = {
                    "status": "match" if match else "differ",
                    "details": [] if match else ["Binary files differ"],
                }

        return results

    def compare_tsv_files(
        self,
        old_file: Path,
        new_file: Path,
        ignore_columns: Optional[List[str]] = None,
        sort_by: Optional[List[str]] = None,
    ) -> Tuple[bool, List[str]]:
        """Compare two TSV files in detail.

        Returns (match, list of differences).
        """
        ignore_columns = ignore_columns or []
        differences = []

        try:
            # Read files
            old_df = pd.read_csv(old_file, sep="\t", dtype=str, na_values=[])
            new_df = pd.read_csv(new_file, sep="\t", dtype=str, na_values=[])

            # Basic shape check
            if old_df.shape[0] != new_df.shape[0]:
                differences.append(
                    f"Row count differs: old={old_df.shape[0]}, new={new_df.shape[0]}"
                )

            if old_df.shape[1] != new_df.shape[1]:
                differences.append(
                    f"Column count differs: old={old_df.shape[1]}, new={new_df.shape[1]}"
                )

            # Column differences
            old_cols = set(old_df.columns) - set(ignore_columns)
            new_cols = set(new_df.columns) - set(ignore_columns)

            missing_cols = old_cols - new_cols
            extra_cols = new_cols - old_cols

            if missing_cols:
                differences.append(f"Missing columns: {sorted(missing_cols)}")
            if extra_cols:
                differences.append(f"Extra columns: {sorted(extra_cols)}")

            # If structure differs significantly, return early
            if missing_cols or extra_cols or old_df.shape[0] != new_df.shape[0]:
                return False, differences

            # Sort for consistent comparison
            if sort_by:
                sort_cols = [c for c in sort_by if c in old_df.columns]
            else:
                # Default sorting by position columns if available
                sort_cols = [c for c in ["CHROM", "POS", "REF", "ALT"] if c in old_df.columns]

            if sort_cols:
                old_df = old_df.sort_values(sort_cols).reset_index(drop=True)
                new_df = new_df.sort_values(sort_cols).reset_index(drop=True)

            # Detailed column comparison
            common_cols = sorted(old_cols & new_cols)
            col_diffs = {}

            for col in common_cols:
                old_vals = old_df[col].fillna("")
                new_vals = new_df[col].fillna("")

                if not old_vals.equals(new_vals):
                    # Count differences
                    diff_mask = old_vals != new_vals
                    n_diff = diff_mask.sum()

                    # Sample differences
                    diff_indices = diff_mask[diff_mask].index[:5]  # First 5
                    samples = []

                    for idx in diff_indices:
                        samples.append(
                            {
                                "row": idx,
                                "old": old_vals.iloc[idx],
                                "new": new_vals.iloc[idx],
                            }
                        )

                    col_diffs[col] = {
                        "count": n_diff,
                        "samples": samples,
                    }

            # Format column differences
            for col, info in col_diffs.items():
                diff_str = f"Column '{col}': {info['count']} differences"
                if self.verbose and info["samples"]:
                    diff_str += "\n  Sample differences:"
                    for s in info["samples"]:
                        diff_str += f"\n    Row {s['row']}: '{s['old']}' → '{s['new']}'"
                differences.append(diff_str)

            return len(col_diffs) == 0, differences

        except Exception as e:
            return False, [f"Error comparing files: {e}"]

    def compare_json_files(
        self,
        old_file: Path,
        new_file: Path,
        ignore_keys: Optional[List[str]] = None,
    ) -> Tuple[bool, List[str]]:
        """Compare two JSON files.

        Returns (match, list of differences).
        """
        ignore_keys = ignore_keys or []
        differences = []

        try:
            with open(old_file) as f:
                old_data = json.load(f)
            with open(new_file) as f:
                new_data = json.load(f)

            # Compare recursively
            self._compare_json_objects(old_data, new_data, "", differences, ignore_keys)

            return len(differences) == 0, differences

        except Exception as e:
            return False, [f"Error comparing JSON files: {e}"]

    def _compare_json_objects(
        self,
        old_obj,
        new_obj,
        path: str,
        differences: List[str],
        ignore_keys: List[str],
    ):
        """Recursively compare JSON objects."""
        if type(old_obj) != type(new_obj):
            differences.append(
                f"{path}: Type mismatch - old={type(old_obj).__name__}, "
                f"new={type(new_obj).__name__}"
            )
            return

        if isinstance(old_obj, dict):
            old_keys = set(old_obj.keys()) - set(ignore_keys)
            new_keys = set(new_obj.keys()) - set(ignore_keys)

            for key in old_keys - new_keys:
                differences.append(f"{path}.{key}: Missing in new")

            for key in new_keys - old_keys:
                differences.append(f"{path}.{key}: Added in new")

            for key in old_keys & new_keys:
                self._compare_json_objects(
                    old_obj[key],
                    new_obj[key],
                    f"{path}.{key}" if path else key,
                    differences,
                    ignore_keys,
                )

        elif isinstance(old_obj, list):
            if len(old_obj) != len(new_obj):
                differences.append(
                    f"{path}: List length differs - old={len(old_obj)}, " f"new={len(new_obj)}"
                )
            else:
                for i, (old_item, new_item) in enumerate(zip(old_obj, new_obj)):
                    self._compare_json_objects(
                        old_item,
                        new_item,
                        f"{path}[{i}]",
                        differences,
                        ignore_keys,
                    )

        elif isinstance(old_obj, (int, float)):
            # Numeric comparison with tolerance
            if abs(old_obj - new_obj) > 1e-6:
                differences.append(f"{path}: Value differs - old={old_obj}, new={new_obj}")

        elif old_obj != new_obj:
            differences.append(f"{path}: Value differs - old={old_obj}, new={new_obj}")

    def _files_identical(self, file1: Path, file2: Path) -> bool:
        """Check if two files are byte-for-byte identical."""
        if file1.stat().st_size != file2.stat().st_size:
            return False

        with open(file1, "rb") as f1, open(file2, "rb") as f2:
            while True:
                b1 = f1.read(8192)
                b2 = f2.read(8192)
                if b1 != b2:
                    return False
                if not b1:
                    return True

    def print_summary(self, results: Dict[str, Dict]):
        """Print a colored summary of comparison results."""
        total = len(results)
        matches = sum(1 for r in results.values() if r["status"] == "match")
        differs = sum(1 for r in results.values() if r["status"] == "differ")
        missing_old = sum(1 for r in results.values() if r["status"] == "missing_in_old")
        missing_new = sum(1 for r in results.values() if r["status"] == "missing_in_new")

        print("\n" + "=" * 60)
        print("COMPARISON SUMMARY")
        print("=" * 60)

        print(f"Total files compared: {total}")
        print(f"{Fore.GREEN}✓ Matching: {matches}{Style.RESET_ALL}")

        if differs > 0:
            print(f"{Fore.RED}✗ Different: {differs}{Style.RESET_ALL}")

        if missing_old > 0:
            print(f"{Fore.YELLOW}+ Only in new: {missing_old}{Style.RESET_ALL}")

        if missing_new > 0:
            print(f"{Fore.YELLOW}- Only in old: {missing_new}{Style.RESET_ALL}")

        print("=" * 60)

        # Detailed results
        if differs > 0 or missing_old > 0 or missing_new > 0:
            print("\nDETAILED RESULTS:")
            print("-" * 60)

            for filename, result in sorted(results.items()):
                if result["status"] == "match":
                    continue

                if result["status"] == "differ":
                    print(f"\n{Fore.RED}✗ {filename}{Style.RESET_ALL}")
                    for detail in result["details"][:10]:  # Limit output
                        print(f"  {detail}")
                    if len(result["details"]) > 10:
                        print(f"  ... and {len(result['details']) - 10} more differences")

                elif result["status"] == "missing_in_old":
                    print(f"\n{Fore.YELLOW}+ {filename} (only in new){Style.RESET_ALL}")

                elif result["status"] == "missing_in_new":
                    print(f"\n{Fore.YELLOW}- {filename} (only in old){Style.RESET_ALL}")


def main():
    """Main entry point."""
    parser = argparse.ArgumentParser(
        description="Compare outputs between old and new pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "old_dir",
        type=Path,
        help="Directory containing old pipeline outputs",
    )

    parser.add_argument(
        "new_dir",
        type=Path,
        help="Directory containing new pipeline outputs",
    )

    parser.add_argument(
        "--pattern",
        default="*.tsv",
        help="File pattern to compare (default: *.tsv)",
    )

    parser.add_argument(
        "--ignore-columns",
        nargs="+",
        help="Columns to ignore in TSV comparison",
    )

    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Show detailed differences",
    )

    parser.add_argument(
        "--output-json",
        type=Path,
        help="Save detailed results to JSON file",
    )

    args = parser.parse_args()

    # Validate directories
    if not args.old_dir.exists():
        print(f"{Fore.RED}Error: Old directory not found: {args.old_dir}{Style.RESET_ALL}")
        sys.exit(1)

    if not args.new_dir.exists():
        print(f"{Fore.RED}Error: New directory not found: {args.new_dir}{Style.RESET_ALL}")
        sys.exit(1)

    # Run comparison
    comparator = OutputComparator(verbose=args.verbose)
    results = comparator.compare_directories(
        args.old_dir,
        args.new_dir,
        args.pattern,
    )

    # Print summary
    comparator.print_summary(results)

    # Save detailed results if requested
    if args.output_json:
        with open(args.output_json, "w") as f:
            json.dump(results, f, indent=2, default=str)
        print(f"\nDetailed results saved to: {args.output_json}")

    # Exit with error if differences found
    if any(r["status"] != "match" for r in results.values()):
        sys.exit(1)


if __name__ == "__main__":
    main()
