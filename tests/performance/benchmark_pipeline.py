#!/usr/bin/env python
"""Performance benchmarking framework for comparing old and new pipeline implementations.

This script measures and compares performance metrics including:
- Execution time
- Memory usage
- CPU utilization
- Scalability with different workloads
"""

import argparse
import json
import os
import platform
import psutil
import subprocess
import sys
import time
from dataclasses import dataclass, asdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


@dataclass
class BenchmarkResult:
    """Container for benchmark results."""

    test_name: str
    pipeline_type: str  # "old" or "new"
    execution_time: float  # seconds
    peak_memory_mb: float
    cpu_percent: float
    n_variants_processed: int
    n_genes: int
    threads: int
    success: bool
    error: Optional[str] = None

    @property
    def variants_per_second(self) -> float:
        """Calculate variant processing rate."""
        if self.execution_time > 0:
            return self.n_variants_processed / self.execution_time
        return 0.0


@dataclass
class BenchmarkConfig:
    """Configuration for a benchmark test."""

    name: str
    gene_file: Optional[Path] = None
    gene_names: Optional[List[str]] = None
    n_genes: int = 1
    preset: Optional[str] = None
    extra_args: Optional[List[str]] = None
    threads: int = 1
    description: str = ""


class PerformanceBenchmark:
    """Runs performance benchmarks on pipelines."""

    def __init__(self, vcf_file: str, output_dir: Path):
        self.vcf_file = vcf_file
        self.output_dir = output_dir
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Create subdirectories
        self.logs_dir = self.output_dir / "logs"
        self.results_dir = self.output_dir / "results"
        self.plots_dir = self.output_dir / "plots"

        for d in [self.logs_dir, self.results_dir, self.plots_dir]:
            d.mkdir(exist_ok=True)

    def run_benchmark(
        self,
        config: BenchmarkConfig,
        use_new_pipeline: bool = False,
    ) -> BenchmarkResult:
        """Run a single benchmark test."""
        pipeline_type = "new" if use_new_pipeline else "old"
        print(f"Running {pipeline_type} pipeline: {config.name}")

        # Prepare output
        output_file = self.results_dir / f"{config.name}_{pipeline_type}.tsv"
        log_file = self.logs_dir / f"{config.name}_{pipeline_type}.log"

        # Build command
        cmd = [
            "variantcentrifuge",
            "--vcf-file",
            self.vcf_file,
            "--output-file",
            str(output_file),
            "--threads",
            str(config.threads),
            "--log-file",
            str(log_file),
        ]

        if use_new_pipeline:
            cmd.append("--use-new-pipeline")

        # Add genes
        if config.gene_file:
            cmd.extend(["--gene-file", str(config.gene_file)])
        elif config.gene_names:
            cmd.extend(["--gene-name", ",".join(config.gene_names)])

        if config.preset:
            cmd.extend(["--preset", config.preset])

        if config.extra_args:
            cmd.extend(config.extra_args)

        # Monitor process
        start_time = time.time()
        peak_memory = 0
        cpu_samples = []

        try:
            # Start process
            proc = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
            )

            # Monitor while running
            psutil_proc = psutil.Process(proc.pid)

            while proc.poll() is None:
                try:
                    # Memory usage (including children)
                    mem_info = psutil_proc.memory_info()
                    current_memory = mem_info.rss / 1024 / 1024  # MB

                    # Include child processes
                    for child in psutil_proc.children(recursive=True):
                        try:
                            child_mem = child.memory_info().rss / 1024 / 1024
                            current_memory += child_mem
                        except (psutil.NoSuchProcess, psutil.AccessDenied):
                            pass

                    peak_memory = max(peak_memory, current_memory)

                    # CPU usage
                    cpu_percent = psutil_proc.cpu_percent(interval=0.1)
                    cpu_samples.append(cpu_percent)

                except (psutil.NoSuchProcess, psutil.AccessDenied):
                    break

                time.sleep(0.1)

            # Get final status
            stdout, stderr = proc.communicate()
            execution_time = time.time() - start_time

            if proc.returncode != 0:
                error_msg = f"Return code: {proc.returncode}\n{stderr}"
                return BenchmarkResult(
                    test_name=config.name,
                    pipeline_type=pipeline_type,
                    execution_time=execution_time,
                    peak_memory_mb=peak_memory,
                    cpu_percent=0.0,
                    n_variants_processed=0,
                    n_genes=config.n_genes,
                    threads=config.threads,
                    success=False,
                    error=error_msg,
                )

            # Count variants processed
            n_variants = self._count_output_variants(output_file)

            # Average CPU usage
            avg_cpu = np.mean(cpu_samples) if cpu_samples else 0.0

            return BenchmarkResult(
                test_name=config.name,
                pipeline_type=pipeline_type,
                execution_time=execution_time,
                peak_memory_mb=peak_memory,
                cpu_percent=avg_cpu,
                n_variants_processed=n_variants,
                n_genes=config.n_genes,
                threads=config.threads,
                success=True,
            )

        except Exception as e:
            return BenchmarkResult(
                test_name=config.name,
                pipeline_type=pipeline_type,
                execution_time=0.0,
                peak_memory_mb=0.0,
                cpu_percent=0.0,
                n_variants_processed=0,
                n_genes=config.n_genes,
                threads=config.threads,
                success=False,
                error=str(e),
            )

    def _count_output_variants(self, output_file: Path) -> int:
        """Count variants in output file."""
        if not output_file.exists():
            return 0

        try:
            # Count lines minus header
            with open(output_file) as f:
                return sum(1 for _ in f) - 1
        except Exception:
            return 0

    def run_comparison(
        self,
        configs: List[BenchmarkConfig],
        runs: int = 3,
    ) -> pd.DataFrame:
        """Run benchmarks for both pipelines and compare.

        Runs each test multiple times and reports statistics.
        """
        results = []

        for config in configs:
            print(f"\nBenchmarking: {config.name}")
            print(f"Description: {config.description}")
            print(f"Genes: {config.n_genes}, Threads: {config.threads}")

            # Run multiple times for each pipeline
            for pipeline in ["old", "new"]:
                use_new = pipeline == "new"

                for run_idx in range(runs):
                    print(f"  {pipeline} pipeline, run {run_idx + 1}/{runs}")

                    result = self.run_benchmark(config, use_new)
                    results.append(
                        {
                            **asdict(result),
                            "run_idx": run_idx,
                        }
                    )

                    if not result.success:
                        print(f"    Failed: {result.error}")
                    else:
                        print(
                            f"    Time: {result.execution_time:.2f}s, "
                            f"Memory: {result.peak_memory_mb:.1f}MB, "
                            f"Variants: {result.n_variants_processed}"
                        )

        # Convert to DataFrame
        df = pd.DataFrame(results)

        # Save raw results
        results_file = self.output_dir / "benchmark_results.csv"
        df.to_csv(results_file, index=False)
        print(f"\nResults saved to: {results_file}")

        return df

    def generate_report(self, df: pd.DataFrame) -> Dict:
        """Generate performance comparison report."""
        report = {
            "timestamp": datetime.now().isoformat(),
            "system_info": {
                "platform": platform.platform(),
                "processor": platform.processor(),
                "cpu_count": os.cpu_count(),
                "memory_gb": psutil.virtual_memory().total / 1024**3,
            },
            "summary": {},
            "details": {},
        }

        # Calculate summary statistics
        summary_stats = []

        for test_name in df["test_name"].unique():
            test_df = df[df["test_name"] == test_name]

            for pipeline in ["old", "new"]:
                pipeline_df = test_df[test_df["pipeline_type"] == pipeline]

                if len(pipeline_df) == 0:
                    continue

                # Only successful runs
                success_df = pipeline_df[pipeline_df["success"]]

                if len(success_df) > 0:
                    stats = {
                        "test_name": test_name,
                        "pipeline": pipeline,
                        "n_runs": len(success_df),
                        "avg_time": success_df["execution_time"].mean(),
                        "std_time": success_df["execution_time"].std(),
                        "min_time": success_df["execution_time"].min(),
                        "max_time": success_df["execution_time"].max(),
                        "avg_memory_mb": success_df["peak_memory_mb"].mean(),
                        "avg_cpu_percent": success_df["cpu_percent"].mean(),
                        "variants_per_sec": success_df["variants_per_second"].mean(),
                    }
                    summary_stats.append(stats)

        summary_df = pd.DataFrame(summary_stats)

        # Calculate speedup
        speedup_results = []
        for test_name in summary_df["test_name"].unique():
            old_stats = summary_df[
                (summary_df["test_name"] == test_name) & (summary_df["pipeline"] == "old")
            ]
            new_stats = summary_df[
                (summary_df["test_name"] == test_name) & (summary_df["pipeline"] == "new")
            ]

            if len(old_stats) > 0 and len(new_stats) > 0:
                speedup = old_stats.iloc[0]["avg_time"] / new_stats.iloc[0]["avg_time"]
                memory_ratio = (
                    new_stats.iloc[0]["avg_memory_mb"] / old_stats.iloc[0]["avg_memory_mb"]
                )

                speedup_results.append(
                    {
                        "test_name": test_name,
                        "speedup": speedup,
                        "memory_ratio": memory_ratio,
                        "old_time": old_stats.iloc[0]["avg_time"],
                        "new_time": new_stats.iloc[0]["avg_time"],
                        "improvement_pct": (1 - 1 / speedup) * 100 if speedup > 0 else 0,
                    }
                )

        report["summary"] = {
            "total_tests": len(df["test_name"].unique()),
            "total_runs": len(df),
            "overall_speedup": np.mean([r["speedup"] for r in speedup_results]),
            "speedup_by_test": speedup_results,
        }

        report["details"] = summary_df.to_dict("records")

        # Save report
        report_file = self.output_dir / "performance_report.json"
        with open(report_file, "w") as f:
            json.dump(report, f, indent=2, default=str)

        return report

    def plot_results(self, df: pd.DataFrame):
        """Generate performance visualization plots."""
        # Setup plot style
        try:
            plt.style.use("seaborn-v0_8-darkgrid")
        except Exception:
            # Fallback to a standard style if seaborn not available
            plt.style.use("ggplot")

        # 1. Execution time comparison
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(15, 6))

        # Bar plot of average execution times
        summary = (
            df.groupby(["test_name", "pipeline_type"])["execution_time"]
            .agg(["mean", "std"])
            .reset_index()
        )

        test_names = summary["test_name"].unique()
        x = np.arange(len(test_names))
        width = 0.35

        old_means = summary[summary["pipeline_type"] == "old"]["mean"].values
        old_stds = summary[summary["pipeline_type"] == "old"]["std"].values
        new_means = summary[summary["pipeline_type"] == "new"]["mean"].values
        new_stds = summary[summary["pipeline_type"] == "new"]["std"].values

        ax1.bar(x - width / 2, old_means, width, yerr=old_stds, label="Old Pipeline", alpha=0.8)
        ax1.bar(x + width / 2, new_means, width, yerr=new_stds, label="New Pipeline", alpha=0.8)

        ax1.set_xlabel("Test Configuration")
        ax1.set_ylabel("Execution Time (seconds)")
        ax1.set_title("Pipeline Execution Time Comparison")
        ax1.set_xticks(x)
        ax1.set_xticklabels(test_names, rotation=45, ha="right")
        ax1.legend()

        # Memory usage comparison
        memory_summary = (
            df.groupby(["test_name", "pipeline_type"])["peak_memory_mb"].mean().reset_index()
        )

        old_memory = memory_summary[memory_summary["pipeline_type"] == "old"][
            "peak_memory_mb"
        ].values
        new_memory = memory_summary[memory_summary["pipeline_type"] == "new"][
            "peak_memory_mb"
        ].values

        ax2.bar(x - width / 2, old_memory, width, label="Old Pipeline", alpha=0.8)
        ax2.bar(x + width / 2, new_memory, width, label="New Pipeline", alpha=0.8)

        ax2.set_xlabel("Test Configuration")
        ax2.set_ylabel("Peak Memory (MB)")
        ax2.set_title("Pipeline Memory Usage Comparison")
        ax2.set_xticks(x)
        ax2.set_xticklabels(test_names, rotation=45, ha="right")
        ax2.legend()

        plt.tight_layout()
        plt.savefig(self.plots_dir / "performance_comparison.png", dpi=300)
        plt.close()

        # 2. Scalability plot (if we have different thread counts)
        thread_tests = df[df["threads"] > 1]
        if len(thread_tests) > 0:
            fig, ax = plt.subplots(figsize=(10, 6))

            for pipeline in ["old", "new"]:
                pipeline_df = thread_tests[thread_tests["pipeline_type"] == pipeline]
                thread_summary = pipeline_df.groupby("threads")["execution_time"].mean()

                ax.plot(
                    thread_summary.index,
                    thread_summary.values,
                    marker="o",
                    linewidth=2,
                    markersize=8,
                    label=f"{pipeline.capitalize()} Pipeline",
                )

            ax.set_xlabel("Number of Threads")
            ax.set_ylabel("Execution Time (seconds)")
            ax.set_title("Multi-threading Performance Scaling")
            ax.legend()
            ax.grid(True, alpha=0.3)

            plt.tight_layout()
            plt.savefig(self.plots_dir / "scalability_plot.png", dpi=300)
            plt.close()

        # 3. Speedup visualization
        speedup_data = []
        for test_name in df["test_name"].unique():
            old_time = df[(df["test_name"] == test_name) & (df["pipeline_type"] == "old")][
                "execution_time"
            ].mean()
            new_time = df[(df["test_name"] == test_name) & (df["pipeline_type"] == "new")][
                "execution_time"
            ].mean()

            if old_time > 0 and new_time > 0:
                speedup = old_time / new_time
                speedup_data.append(
                    {
                        "test": test_name,
                        "speedup": speedup,
                    }
                )

        if speedup_data:
            speedup_df = pd.DataFrame(speedup_data)

            fig, ax = plt.subplots(figsize=(10, 6))
            bars = ax.bar(speedup_df["test"], speedup_df["speedup"])

            # Color bars based on speedup
            for i, (bar, speedup) in enumerate(zip(bars, speedup_df["speedup"])):
                if speedup > 1.2:
                    bar.set_color("green")
                elif speedup > 0.9:
                    bar.set_color("orange")
                else:
                    bar.set_color("red")

            ax.axhline(y=1.0, color="black", linestyle="--", alpha=0.5)
            ax.set_xlabel("Test Configuration")
            ax.set_ylabel("Speedup Factor")
            ax.set_title("New Pipeline Speedup vs Old Pipeline")
            ax.set_xticklabels(speedup_df["test"], rotation=45, ha="right")

            # Add value labels on bars
            for bar, speedup in zip(bars, speedup_df["speedup"]):
                height = bar.get_height()
                ax.text(
                    bar.get_x() + bar.get_width() / 2.0,
                    height,
                    f"{speedup:.2f}x",
                    ha="center",
                    va="bottom",
                )

            plt.tight_layout()
            plt.savefig(self.plots_dir / "speedup_comparison.png", dpi=300)
            plt.close()


def create_benchmark_configs() -> List[BenchmarkConfig]:
    """Create standard benchmark configurations."""
    configs = [
        # Single gene tests
        BenchmarkConfig(
            name="single_gene_1thread",
            gene_names=["BRCA1"],
            n_genes=1,
            threads=1,
            description="Single gene, single thread baseline",
        ),
        BenchmarkConfig(
            name="single_gene_4thread",
            gene_names=["BRCA1"],
            n_genes=1,
            threads=4,
            description="Single gene with parallelization",
        ),
        # Small gene set
        BenchmarkConfig(
            name="small_geneset_1thread",
            gene_names=["BRCA1", "BRCA2", "TP53", "EGFR", "KRAS"],
            n_genes=5,
            threads=1,
            description="Small gene set (5 genes)",
        ),
        BenchmarkConfig(
            name="small_geneset_4thread",
            gene_names=["BRCA1", "BRCA2", "TP53", "EGFR", "KRAS"],
            n_genes=5,
            threads=4,
            description="Small gene set with parallelization",
        ),
        # Medium gene set
        BenchmarkConfig(
            name="medium_geneset_1thread",
            n_genes=50,
            threads=1,
            description="Medium gene set (50 genes)",
        ),
        BenchmarkConfig(
            name="medium_geneset_8thread",
            n_genes=50,
            threads=8,
            description="Medium gene set with high parallelization",
        ),
        # Complex filtering
        BenchmarkConfig(
            name="complex_filter",
            gene_names=["BRCA1"],
            n_genes=1,
            preset="rare,coding,pathogenic",
            threads=1,
            description="Complex filtering pipeline",
        ),
        # With scoring
        BenchmarkConfig(
            name="with_scoring",
            gene_names=["CFTR"],
            n_genes=1,
            extra_args=["--scoring-config-path", "scoring/inheritance_score"],
            threads=1,
            description="Pipeline with variant scoring",
        ),
    ]

    return configs


def main():
    """Run the benchmark pipeline performance comparison."""
    parser = argparse.ArgumentParser(
        description="Benchmark pipeline performance",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "vcf_file",
        help="VCF file to use for benchmarking",
    )

    parser.add_argument(
        "--output-dir",
        default="benchmark_results",
        help="Output directory for results (default: benchmark_results)",
    )

    parser.add_argument(
        "--runs",
        type=int,
        default=3,
        help="Number of runs per test (default: 3)",
    )

    parser.add_argument(
        "--config",
        choices=["quick", "standard", "extensive"],
        default="standard",
        help="Benchmark configuration set (default: standard)",
    )

    args = parser.parse_args()

    # Check VCF exists
    if not Path(args.vcf_file).exists():
        print(f"Error: VCF file not found: {args.vcf_file}")
        sys.exit(1)

    # Select configurations
    all_configs = create_benchmark_configs()

    if args.config == "quick":
        # Just single gene tests
        configs = [c for c in all_configs if "single_gene" in c.name]
    elif args.config == "extensive":
        # All tests plus additional large ones
        configs = all_configs + [
            BenchmarkConfig(
                name="large_geneset",
                n_genes=200,
                threads=8,
                description="Large gene set (200 genes)",
            ),
        ]
    else:
        # Standard set
        configs = all_configs

    # Prepare gene files for tests that need them
    output_dir = Path(args.output_dir)
    gene_files_dir = output_dir / "gene_files"
    gene_files_dir.mkdir(parents=True, exist_ok=True)

    # Create gene files for tests that specify n_genes but no gene names
    cancer_genes = [
        "ABL1",
        "AKT1",
        "ALK",
        "APC",
        "AR",
        "ATM",
        "BRAF",
        "BRCA1",
        "BRCA2",
        "CDH1",
        "CDKN2A",
        "CTNNB1",
        "EGFR",
        "ERBB2",
        "EZH2",
        "FLT3",
        "GATA3",
        "HER2",
        "IDH1",
        "JAK2",
        "KIT",
        "KRAS",
        "MET",
        "MLH1",
        "MYC",
        "NF1",
        "NOTCH1",
        "NPM1",
        "NRAS",
        "PALB2",
        "PDGFRA",
        "PIK3CA",
        "PTEN",
        "PTPN11",
        "RB1",
        "RET",
        "RUNX1",
        "SMAD4",
        "STK11",
        "TP53",
        "TSC1",
        "VHL",
        "WT1",
    ]

    for config in configs:
        if config.n_genes > 1 and not config.gene_names and not config.gene_file:
            # Create gene file
            gene_file = gene_files_dir / f"{config.name}_genes.txt"
            # Repeat genes if needed
            genes_to_use = []
            while len(genes_to_use) < config.n_genes:
                genes_to_use.extend(cancer_genes)
            genes_to_use = genes_to_use[: config.n_genes]

            gene_file.write_text("\n".join(genes_to_use))
            config.gene_file = gene_file

    # Run benchmarks
    print(f"Running {len(configs)} benchmark configurations")
    print(f"Each test will run {args.runs} times")
    print(f"Output directory: {output_dir}")
    print()

    benchmark = PerformanceBenchmark(args.vcf_file, output_dir)

    # Run comparison
    results_df = benchmark.run_comparison(configs, runs=args.runs)

    # Generate report
    report = benchmark.generate_report(results_df)

    # Generate plots
    benchmark.plot_results(results_df)

    # Print summary
    print("\n" + "=" * 60)
    print("BENCHMARK SUMMARY")
    print("=" * 60)

    print(f"\nOverall speedup: {report['summary']['overall_speedup']:.2f}x")

    print("\nSpeedup by test:")
    for test_result in report["summary"]["speedup_by_test"]:
        print(
            f"  {test_result['test_name']}: {test_result['speedup']:.2f}x "
            f"({test_result['improvement_pct']:.1f}% improvement)"
        )

    print(f"\nDetailed results saved to: {output_dir}")
    print("  - Raw data: benchmark_results.csv")
    print("  - Report: performance_report.json")
    print("  - Plots: plots/")


if __name__ == "__main__":
    main()
