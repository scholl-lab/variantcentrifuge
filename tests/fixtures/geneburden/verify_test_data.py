#!/usr/bin/env python3
"""
Verification script for comprehensive gene burden test data.

This script verifies that the generated test data is properly formatted
and ready for gene burden analysis testing.
"""

import json
import sys
from pathlib import Path
import gzip

def verify_test_data(data_dir: Path) -> bool:
    """Verify all components of the test dataset."""
    print(f"Verifying test dataset in: {data_dir}")
    
    success = True
    
    # Check core files
    required_files = [
        "test_data.vcf.gz",
        "case_samples.txt", 
        "control_samples.txt",
        "test_genes.txt",
        "phenotypes_basic.csv",
        "dataset_statistics.json",
        "README.md"
    ]
    
    print("\n1. Checking required files...")
    for filename in required_files:
        filepath = data_dir / filename
        if filepath.exists():
            print(f"  ✓ {filename}")
        else:
            print(f"  ✗ {filename} - MISSING")
            success = False
    
    # Check VCF structure
    print("\n2. Checking VCF structure...")
    vcf_path = data_dir / "test_data.vcf.gz"
    if vcf_path.exists():
        try:
            with gzip.open(vcf_path, 'rt') as f:
                header_found = False
                variant_count = 0
                sample_count = 0
                
                for line in f:
                    if line.startswith('##'):
                        continue
                    elif line.startswith('#CHROM'):
                        header_found = True
                        fields = line.strip().split('\t')
                        sample_count = len(fields) - 9  # Subtract standard VCF columns
                        print(f"  ✓ VCF header found with {sample_count} samples")
                    else:
                        variant_count += 1
                        if variant_count > 5:  # Just check first few variants
                            break
                
                if header_found and variant_count > 0:
                    print(f"  ✓ VCF contains {variant_count}+ variants")
                else:
                    print("  ✗ VCF structure invalid")
                    success = False
                    
        except Exception as e:
            print(f"  ✗ Error reading VCF: {e}")
            success = False
    
    # Check sample files
    print("\n3. Checking sample files...")
    case_file = data_dir / "case_samples.txt"
    control_file = data_dir / "control_samples.txt"
    
    if case_file.exists() and control_file.exists():
        case_samples = case_file.read_text().strip().split('\n')
        control_samples = control_file.read_text().strip().split('\n')
        
        print(f"  ✓ Case samples: {len(case_samples)}")
        print(f"  ✓ Control samples: {len(control_samples)}")
        
        # Check for overlap
        case_set = set(case_samples)
        control_set = set(control_samples)
        overlap = case_set & control_set
        
        if overlap:
            print(f"  ⚠ Sample overlap detected: {overlap}")
        else:
            print("  ✓ No sample overlap between cases and controls")
    
    # Check phenotype files
    print("\n4. Checking phenotype files...")
    pheno_files = [
        "phenotypes_basic.csv",
        "phenotypes_extended.csv", 
        "phenotypes.tsv",
        "phenotypes_alt_columns.csv"
    ]
    
    for pheno_file in pheno_files:
        filepath = data_dir / pheno_file
        if filepath.exists():
            try:
                with open(filepath, 'r') as f:
                    lines = f.readlines()
                    print(f"  ✓ {pheno_file}: {len(lines)-1} data rows")
            except Exception as e:
                print(f"  ✗ Error reading {pheno_file}: {e}")
                success = False
    
    # Check gene files
    print("\n5. Checking gene files...")
    gene_files = ["test_genes.txt", "disease_genes.txt", "control_genes.txt"]
    
    for gene_file in gene_files:
        filepath = data_dir / gene_file
        if filepath.exists():
            genes = filepath.read_text().strip().split('\n')
            print(f"  ✓ {gene_file}: {len(genes)} genes")
    
    # Check HPO term files
    print("\n6. Checking HPO term files...")
    hpo_files = ["case_hpo_terms.txt", "control_hpo_terms.txt", "mixed_hpo_terms.txt"]
    
    for hpo_file in hpo_files:
        filepath = data_dir / hpo_file
        if filepath.exists():
            terms = filepath.read_text().strip().split('\n')
            print(f"  ✓ {hpo_file}: {len(terms)} HPO terms")
    
    # Check statistics
    print("\n7. Checking dataset statistics...")
    stats_file = data_dir / "dataset_statistics.json"
    if stats_file.exists():
        try:
            with open(stats_file, 'r') as f:
                stats = json.load(f)
                
            print(f"  ✓ Total samples: {stats['sample_composition']['total_samples']}")
            print(f"  ✓ Case samples: {stats['sample_composition']['case_samples']}")
            print(f"  ✓ Control samples: {stats['sample_composition']['control_samples']}")
            print(f"  ✓ Test genes: {len(stats['test_genes']['disease_genes']) + len(stats['test_genes']['control_genes'])}")
            
        except Exception as e:
            print(f"  ✗ Error reading statistics: {e}")
            success = False
    
    # Final summary
    print("\n" + "="*50)
    if success:
        print("✅ ALL CHECKS PASSED - Test dataset ready for use!")
        print(f"\nNext steps:")
        print(f"1. cd {data_dir}")
        print(f"2. ./run_comprehensive_tests.sh")
    else:
        print("❌ SOME CHECKS FAILED - Please review errors above")
    
    return success


def main():
    """Main verification function."""
    if len(sys.argv) > 1:
        data_dir = Path(sys.argv[1])
    else:
        # Default to the generated test data directory
        script_dir = Path(__file__).parent
        data_dir = script_dir / "comprehensive_test_data"
    
    if not data_dir.exists():
        print(f"Error: Test data directory not found: {data_dir}")
        sys.exit(1)
    
    success = verify_test_data(data_dir)
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()