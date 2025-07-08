#!/usr/bin/env python
"""
Demo script showing the performance benefits of parallel inheritance analysis.
"""

import time
import pandas as pd
import numpy as np
from variantcentrifuge.inheritance.analyzer import analyze_inheritance
from variantcentrifuge.inheritance.parallel_analyzer import analyze_inheritance_parallel

# Create synthetic data with multiple genes
def create_test_data(n_variants=1000, n_genes=50, n_samples=10):
    """Create test data with variants across multiple genes."""
    np.random.seed(42)
    
    # Create base columns
    data = {
        "CHROM": [f"chr{np.random.randint(1, 23)}" for _ in range(n_variants)],
        "POS": np.random.randint(1000000, 50000000, n_variants),
        "REF": np.random.choice(["A", "T", "G", "C"], n_variants),
        "ALT": np.random.choice(["A", "T", "G", "C"], n_variants),
        "GENE": [f"GENE{i % n_genes}" for i in range(n_variants)],
    }
    
    # Create genotypes for samples
    for i in range(n_samples):
        sample_id = f"sample_{i}"
        # Create realistic genotype patterns
        genotypes = []
        for j in range(n_variants):
            r = np.random.random()
            if r < 0.7:
                gt = "0/0"  # Reference homozygous
            elif r < 0.9:
                gt = "0/1"  # Heterozygous
            else:
                gt = "1/1"  # Alternate homozygous
            genotypes.append(gt)
        data[sample_id] = genotypes
    
    return pd.DataFrame(data)

# Create pedigree data (trio)
def create_pedigree(n_samples=10):
    """Create pedigree data for testing."""
    pedigree = {}
    
    # Create families
    for i in range(0, n_samples, 3):
        if i + 2 < n_samples:
            # Proband
            pedigree[f"sample_{i}"] = {
                "sample_id": f"sample_{i}",
                "father_id": f"sample_{i+1}",
                "mother_id": f"sample_{i+2}",
                "sex": "1",
                "affected_status": "2",
            }
            # Father
            pedigree[f"sample_{i+1}"] = {
                "sample_id": f"sample_{i+1}",
                "father_id": "0",
                "mother_id": "0",
                "sex": "1",
                "affected_status": "1",
            }
            # Mother
            pedigree[f"sample_{i+2}"] = {
                "sample_id": f"sample_{i+2}",
                "father_id": "0",
                "mother_id": "0",
                "sex": "2",
                "affected_status": "1",
            }
    
    return pedigree

def main():
    """Run performance comparison."""
    print("Parallel Inheritance Analysis Performance Demo")
    print("=" * 60)
    
    # Test with different dataset sizes
    test_sizes = [
        (100, 10, 3),    # Small: 100 variants, 10 genes, 3 samples
        (1000, 50, 6),   # Medium: 1000 variants, 50 genes, 6 samples
        (5000, 100, 9),  # Large: 5000 variants, 100 genes, 9 samples
    ]
    
    for n_variants, n_genes, n_samples in test_sizes:
        print(f"\nDataset: {n_variants} variants, {n_genes} genes, {n_samples} samples")
        print("-" * 50)
        
        # Create test data
        df = create_test_data(n_variants, n_genes, n_samples)
        pedigree = create_pedigree(n_samples)
        sample_list = [f"sample_{i}" for i in range(n_samples)]
        
        # Run sequential analysis
        print("Running sequential analysis...")
        start = time.time()
        result_seq = analyze_inheritance(df, pedigree, sample_list)
        seq_time = time.time() - start
        print(f"Sequential time: {seq_time:.2f}s")
        
        # Run parallel analysis with different worker counts
        for n_workers in [2, 4]:
            print(f"\nRunning parallel analysis with {n_workers} workers...")
            start = time.time()
            result_par = analyze_inheritance_parallel(
                df, pedigree, sample_list,
                n_workers=n_workers,
                min_variants_for_parallel=50
            )
            par_time = time.time() - start
            print(f"Parallel time: {par_time:.2f}s")
            print(f"Speedup: {seq_time/par_time:.2f}x")
            
            # Verify results are the same
            patterns_match = result_seq["Inheritance_Pattern"].equals(
                result_par["Inheritance_Pattern"]
            )
            print(f"Results match: {patterns_match}")
        
        # Show pattern distribution
        pattern_counts = result_seq["Inheritance_Pattern"].value_counts()
        print(f"\nPattern distribution:")
        for pattern, count in pattern_counts.items():
            print(f"  {pattern}: {count}")

if __name__ == "__main__":
    main()