"""
Pytest-benchmark result comparison utility.

Simple tool for comparing two benchmark result JSON files and showing what
got faster/slower. Designed for quick iteration during optimization sprints.

Usage
-----
python -m tests.performance.helpers.result_diff baseline.json current.json
"""

import json
import sys
from pathlib import Path


# ANSI color codes
GREEN = "\033[92m"
RED = "\033[91m"
YELLOW = "\033[93m"
RESET = "\033[0m"


def load_benchmark_json(path: str) -> dict:
    """
    Load and validate a pytest-benchmark JSON result file.

    Parameters
    ----------
    path : str
        Path to JSON file

    Returns
    -------
    dict
        Parsed JSON data with 'benchmarks' key

    Raises
    ------
    FileNotFoundError
        If file doesn't exist
    ValueError
        If JSON structure is invalid (missing 'benchmarks' key)
    """
    file_path = Path(path)

    if not file_path.exists():
        raise FileNotFoundError(
            f"Benchmark file not found: {path}\n"
            f"Make sure you've run benchmarks with --benchmark-save first:\n"
            f"  pytest tests/performance/ --benchmark-save=baseline"
        )

    with open(file_path, "r", encoding="utf-8") as f:
        data = json.load(f)

    if "benchmarks" not in data:
        raise ValueError(
            f"Invalid benchmark JSON file: {path}\n"
            f"Expected 'benchmarks' key but found: {list(data.keys())}"
        )

    return data


def compare_benchmark_results(baseline_path: str, current_path: str) -> None:
    """
    Compare two pytest-benchmark JSON files and print a comparison table.

    Shows benchmark name, baseline time, current time, percentage change,
    and status (FASTER/SLOWER/SAME). Uses color coding for easy visual parsing.

    Parameters
    ----------
    baseline_path : str
        Path to baseline benchmark JSON file
    current_path : str
        Path to current benchmark JSON file
    """
    # Load both files
    baseline_data = load_benchmark_json(baseline_path)
    current_data = load_benchmark_json(current_path)

    # Extract benchmark results
    baseline_benchmarks = {b["name"]: b for b in baseline_data["benchmarks"]}
    current_benchmarks = {b["name"]: b for b in current_data["benchmarks"]}

    # Find matched, new, and removed benchmarks
    matched_names = set(baseline_benchmarks.keys()) & set(current_benchmarks.keys())
    new_names = set(current_benchmarks.keys()) - set(baseline_benchmarks.keys())
    removed_names = set(baseline_benchmarks.keys()) - set(current_benchmarks.keys())

    # Print header
    print("\n" + "=" * 100)
    print("BENCHMARK COMPARISON")
    print("=" * 100)
    print(
        f"{'Benchmark':<50} {'Baseline':>12} {'Current':>12} {'Change':>10} {'Status':>10}"
    )
    print("-" * 100)

    # Track statistics
    faster_count = 0
    slower_count = 0
    unchanged_count = 0

    # Compare matched benchmarks
    for name in sorted(matched_names):
        baseline_bench = baseline_benchmarks[name]
        current_bench = current_benchmarks[name]

        # Extract mean times (pytest-benchmark stores this in stats.mean)
        baseline_mean = baseline_bench["stats"]["mean"]
        current_mean = current_bench["stats"]["mean"]

        # Calculate percentage change
        # Positive = slower (bad), Negative = faster (good)
        pct_change = ((current_mean - baseline_mean) / baseline_mean) * 100

        # Determine status
        if abs(pct_change) < 5.0:
            status = "SAME"
            color = RESET
            unchanged_count += 1
        elif pct_change < 0:
            status = f"{GREEN}FASTER{RESET}"
            color = GREEN
            faster_count += 1
        else:
            status = f"{RED}SLOWER{RESET}"
            color = RED
            slower_count += 1

        # Format times (convert seconds to ms for readability)
        baseline_str = f"{baseline_mean * 1000:.2f}ms"
        current_str = f"{current_mean * 1000:.2f}ms"
        change_str = f"{color}{pct_change:+.1f}%{RESET}"

        # Truncate long benchmark names
        display_name = name[:47] + "..." if len(name) > 50 else name

        print(f"{display_name:<50} {baseline_str:>12} {current_str:>12} {change_str:>18} {status}")

    # Print new benchmarks
    if new_names:
        print("\n" + "-" * 100)
        print(f"{YELLOW}NEW BENCHMARKS (in current, not in baseline):{RESET}")
        for name in sorted(new_names):
            bench = current_benchmarks[name]
            mean = bench["stats"]["mean"]
            display_name = name[:47] + "..." if len(name) > 50 else name
            print(f"  {display_name:<50} {mean * 1000:>12.2f}ms")

    # Print removed benchmarks
    if removed_names:
        print("\n" + "-" * 100)
        print(f"{YELLOW}REMOVED BENCHMARKS (in baseline, not in current):{RESET}")
        for name in sorted(removed_names):
            bench = baseline_benchmarks[name]
            mean = bench["stats"]["mean"]
            display_name = name[:47] + "..." if len(name) > 50 else name
            print(f"  {display_name:<50} {mean * 1000:>12.2f}ms")

    # Print summary
    print("\n" + "=" * 100)
    print("SUMMARY")
    print("-" * 100)
    print(
        f"{GREEN}{faster_count} faster{RESET}, "
        f"{RED}{slower_count} slower{RESET}, "
        f"{unchanged_count} unchanged"
    )
    if new_names:
        print(f"{YELLOW}{len(new_names)} new{RESET}")
    if removed_names:
        print(f"{YELLOW}{len(removed_names)} removed{RESET}")
    print("=" * 100 + "\n")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python -m tests.performance.helpers.result_diff baseline.json current.json")
        print("\nExample:")
        print("  # Save baseline")
        print("  pytest tests/performance/ --benchmark-save=baseline")
        print("  # Make changes, save new results")
        print("  pytest tests/performance/ --benchmark-save=current")
        print("  # Compare")
        print("  python -m tests.performance.helpers.result_diff \\")
        print("    .benchmarks/.../baseline.json .benchmarks/.../current.json")
        sys.exit(1)

    baseline_path = sys.argv[1]
    current_path = sys.argv[2]

    try:
        compare_benchmark_results(baseline_path, current_path)
    except (FileNotFoundError, ValueError) as e:
        print(f"\n{RED}Error:{RESET} {e}\n", file=sys.stderr)
        sys.exit(1)
