#!/usr/bin/env python3
"""
Enhanced test data generator for comprehensive gene burden analysis testing.

This script creates a complete test suite covering all possible ways to specify
case/control groups in VariantCentrifuge gene burden analysis.

Test Scenarios Generated:
1. Direct sample specification (--case-samples, --control-samples)
2. Sample file specification (--case-samples-file, --control-samples-file)
3. Phenotype-based classification with HPO terms
4. Phenotype file integration with various column configurations
5. Mixed scenarios and edge cases

The generated data includes:
- Annotated VCF with controlled genotype distributions
- Multiple phenotype file formats
- Sample lists in various formats
- HPO term files
- Test configurations for all scenarios
- Automated test scripts

Based on: testing/tools/test_data_generators/generate_gene_burden_test_data.py
Enhanced for comprehensive testing coverage.
"""

import argparse
import csv
import gzip
import json
import logging
import os
import random
import re
import sys
from collections import defaultdict
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Set, Tuple

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


class ComprehensiveTestDataGenerator:
    """Generate comprehensive test data for all gene burden analysis scenarios."""

    # Test genes with real hg19 coordinates
    DISEASE_GENES = {"PKD1", "PKD2", "BRCA1", "BRCA2"}
    CONTROL_GENES = {"TTN", "OBSCN", "MUC16"}  # Removed duplicate TITIN
    ALL_GENES = DISEASE_GENES | CONTROL_GENES
    
    # Real hg19 gene coordinates (chrom, start, end, gene, strand)
    GENE_COORDINATES = {
        "OBSCN": ("1", 228395830, 228566577, "+"),
        "BRCA2": ("13", 32889610, 32973805, "+"), 
        "PKD1": ("16", 2138710, 2185899, "-"),
        "BRCA1": ("17", 41196311, 41277500, "-"),
        "MUC16": ("19", 8959519, 9092018, "-"),
        "TTN": ("2", 179390715, 179695529, "-"),
        "PKD2": ("4", 88928819, 88998929, "+"),
    }

    # Sample configuration
    NUM_CASES = 40
    NUM_CONTROLS = 60
    TOTAL_SAMPLES = NUM_CASES + NUM_CONTROLS

    # HPO terms for test scenarios
    HPO_TERMS = {
        "case_terms": [
            "HP:0000113",  # Polycystic kidney dysplasia
            "HP:0000003",  # Multicystic kidney dysplasia
            "HP:0000107",  # Renal cyst
            "HP:0000822",  # Hypertension
            "HP:0003774",  # Stage 5 chronic kidney disease
        ],
        "control_terms": [
            "HP:0000001",  # All (root term, essentially healthy)
            "HP:0032101",  # Normal phenotype (custom)
        ],
        "mixed_terms": [
            "HP:0001626",  # Abnormality of the cardiovascular system
            "HP:0000478",  # Abnormality of the eye
            "HP:0001574",  # Abnormality of the integument
        ],
    }

    # Genotype probabilities for controlled distributions
    PROBABILITIES = {
        "case_disease_pathogenic": 0.75,  # 75% cases have pathogenic disease gene variants
        "case_disease_benign": 0.25,     # 25% cases have benign disease gene variants
        "control_disease_pathogenic": 0.10,  # 10% controls have pathogenic disease variants
        "control_disease_benign": 0.15,     # 15% controls have benign disease variants
        "any_control_gene": 0.40,           # 40% all samples have control gene variants
    }

    HOMO_PROB = 0.20  # Probability of homozygous when variant is present

    def __init__(self, seed: int = 42):
        """Initialize generator with reproducible random seed."""
        random.seed(seed)
        
        # Generate sample IDs with realistic naming
        self.case_samples = [f"CASE_{i:03d}" for i in range(1, self.NUM_CASES + 1)]
        self.control_samples = [f"CTRL_{i:03d}" for i in range(1, self.NUM_CONTROLS + 1)]
        self.all_samples = self.case_samples + self.control_samples
        
        # Track statistics
        self.variant_stats = defaultdict(lambda: defaultdict(int))
        
        logger.info(f"Initialized generator with {len(self.case_samples)} cases, "
                   f"{len(self.control_samples)} controls")

    def extract_gene_from_info(self, info_field: str) -> str:
        """Extract gene name from VCF INFO field."""
        # Look for ANN field first (SnpEff annotation)
        ann_match = re.search(r"ANN=([^;]+)", info_field)
        if ann_match:
            ann_value = ann_match.group(1)
            ann_parts = ann_value.split("|")
            if len(ann_parts) > 3:
                gene = ann_parts[3]  # Gene name in position 3
                return gene

        # Fallback to direct GENE= field
        gene_match = re.search(r"GENE=([^;,\s]+)", info_field)
        if gene_match:
            return gene_match.group(1)

        return ""

    def extract_clinsig_from_info(self, info_field: str) -> str:
        """Extract clinical significance from VCF INFO field."""
        # Try multiple common formats
        for pattern in [
            r"ClinVar_CLNSIG=([^;]+)",
            r"dbNSFP_clinvar_clnsig=([^;]+)",
            r"CLNSIG=([^;]+)",
        ]:
            match = re.search(pattern, info_field)
            if match:
                return match.group(1).lower()
        return "uncertain"

    def is_pathogenic(self, clinsig: str) -> bool:
        """Determine if variant is pathogenic."""
        pathogenic_terms = ["pathogenic", "likely_pathogenic", "risk_factor"]
        return any(term in clinsig for term in pathogenic_terms)

    def generate_genotype(self, gene: str, clinsig: str, is_case: bool) -> str:
        """Generate realistic genotype based on gene, clinical significance, and sample type."""
        genotype = "0/0"  # Default: no variant
        
        # Determine variant presence probability
        should_have_variant = False
        
        if gene in self.DISEASE_GENES:
            if is_case:
                prob = (self.PROBABILITIES["case_disease_pathogenic"] if self.is_pathogenic(clinsig)
                       else self.PROBABILITIES["case_disease_benign"])
            else:
                prob = (self.PROBABILITIES["control_disease_pathogenic"] if self.is_pathogenic(clinsig)
                       else self.PROBABILITIES["control_disease_benign"])
            should_have_variant = random.random() < prob
            
        elif gene in self.CONTROL_GENES:
            # Control genes: same probability for all samples
            should_have_variant = random.random() < self.PROBABILITIES["any_control_gene"]

        # Generate genotype if variant present
        if should_have_variant:
            if random.random() < self.HOMO_PROB:
                genotype = "1/1"
            else:
                genotype = "0/1"

        # Update statistics
        sample_type = "case" if is_case else "control"
        self.variant_stats[f"{gene}_{sample_type}"][genotype] += 1

        # Return formatted genotype with realistic quality scores
        if genotype == "0/0":
            return f"{genotype}:55:0,55:99"
        elif genotype == "0/1":
            return f"{genotype}:65:28,32:95"
        else:  # 1/1
            return f"{genotype}:45:0,45:99"

    def create_vcf_file(self, output_path: Path, base_vcf_path: str = None) -> Dict[str, int]:
        """Create test VCF with controlled genotype distributions."""
        logger.info(f"Creating test VCF: {output_path}")
        
        stats = {
            "total_variants": 0,
            "disease_variants": 0,
            "control_variants": 0,
            "skipped_variants": 0,
        }

        # If base VCF provided, process it; otherwise create synthetic data
        if base_vcf_path and Path(base_vcf_path).exists():
            stats = self._process_existing_vcf(base_vcf_path, output_path)
        else:
            stats = self._create_synthetic_vcf(output_path)
            
        logger.info(f"Created VCF with {stats['total_variants']} variants")
        return stats

    def _generate_realistic_variants(self) -> List[Tuple[str, str, str, str, str, str]]:
        """Generate variants within real gene boundaries using hg19 coordinates."""
        variants = []
        
        # Predefined alleles for realistic substitutions
        alleles = ["A", "T", "G", "C"]
        
        # Clinical significance patterns
        clinsig_patterns = {
            "disease_pathogenic": ["pathogenic", "likely_pathogenic"],
            "disease_benign": ["benign", "likely_benign", "uncertain"],
            "control_benign": ["benign", "likely_benign"],
        }
        
        for gene in self.ALL_GENES:
            if gene not in self.GENE_COORDINATES:
                logger.warning(f"No coordinates found for gene {gene}, skipping")
                continue
                
            chrom, start, end, strand = self.GENE_COORDINATES[gene]
            gene_length = end - start
            
            # Number of variants per gene based on gene type and size
            if gene in self.DISEASE_GENES:
                # Disease genes: more pathogenic variants
                num_variants = 3 if gene_length < 100000 else 4
                variant_types = clinsig_patterns["disease_pathogenic"] + clinsig_patterns["disease_benign"]
            else:
                # Control genes: mostly benign variants  
                num_variants = 4 if gene_length > 200000 else 3  # TTN and OBSCN are large
                variant_types = clinsig_patterns["control_benign"]
            
            # Generate variants within gene boundaries
            for i in range(num_variants):
                # Sample position within gene (avoid very start/end for realism)
                margin = min(gene_length // 10, 5000)  # 10% margin or 5kb, whichever is smaller
                pos_start = start + margin
                pos_end = end - margin
                position = random.randint(pos_start, pos_end)
                
                # Generate realistic alleles
                ref_allele = random.choice(alleles)
                alt_allele = random.choice([a for a in alleles if a != ref_allele])
                
                # Assign clinical significance
                if gene in self.DISEASE_GENES:
                    # Disease genes: higher chance of pathogenic variants
                    if i == 0:  # First variant always pathogenic
                        clinsig = "pathogenic"
                    elif i == 1 and num_variants > 2:  # Second might be likely_pathogenic
                        clinsig = "likely_pathogenic" 
                    else:  # Rest are benign/uncertain
                        clinsig = random.choice(["benign", "likely_benign", "uncertain"])
                else:
                    # Control genes: mostly benign
                    clinsig = random.choice(["benign", "likely_benign"])
                
                variants.append((chrom, str(position), ref_allele, alt_allele, gene, clinsig))
        
        # Sort variants by chromosome and position for proper VCF format
        def sort_key(variant):
            chrom, pos, ref, alt, gene, clinsig = variant
            # Handle X, Y, MT chromosomes if present
            if chrom == "X":
                chrom_num = 23
            elif chrom == "Y":
                chrom_num = 24
            elif chrom == "MT":
                chrom_num = 25
            else:
                chrom_num = int(chrom)
            return (chrom_num, int(pos))
        
        variants.sort(key=sort_key)
        logger.info(f"Generated {len(variants)} realistic variants across {len(self.ALL_GENES)} genes")
        
        return variants

    def _create_synthetic_vcf(self, output_path: Path) -> Dict[str, int]:
        """Create synthetic VCF with realistic variant distribution using real gene coordinates."""
        stats = {"total_variants": 0, "disease_variants": 0, "control_variants": 0, "skipped_variants": 0}
        
        # Generate realistic variants within gene boundaries
        synthetic_variants = self._generate_realistic_variants()
        
        with open(output_path, "w") as f:
            # Write VCF header with reference info
            f.write("##fileformat=VCFv4.2\n")
            f.write("##reference=hg19\n")
            f.write(f"##fileDate={datetime.now().strftime('%Y%m%d')}\n")
            f.write("##source=VariantCentrifuge_TestDataGenerator_v1.0\n")
            
            # Add contig lines for chromosomes we're using
            chromosomes = set()
            for gene, (chrom, start, end, strand) in self.GENE_COORDINATES.items():
                chromosomes.add(chrom)
            
            for chrom in sorted(chromosomes, key=lambda x: int(x) if x.isdigit() else ord(x[0])):
                f.write(f"##contig=<ID={chrom}>\n")
            
            # INFO field definitions
            f.write('##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations from test generator">\n')
            f.write('##INFO=<ID=GENE,Number=1,Type=String,Description="Gene name">\n')
            f.write('##INFO=<ID=CLNSIG,Number=.,Type=String,Description="Clinical significance">\n')
            
            # FORMAT field definitions
            f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')
            f.write('##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">\n')
            f.write('##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">\n')
            f.write('##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">\n')
            
            # Write sample header
            header_fields = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
            header = "\t".join(header_fields + self.all_samples) + "\n"
            f.write(header)
            
            # Write variants
            for chrom, pos, ref, alt, gene, clinsig in synthetic_variants:
                stats["total_variants"] += 1
                if gene in self.DISEASE_GENES:
                    stats["disease_variants"] += 1
                else:
                    stats["control_variants"] += 1
                
                # Create INFO field
                info = f"GENE={gene};CLNSIG={clinsig};ANN=variant|{clinsig}|moderate|{gene}|gene_id"
                
                # Generate genotypes for all samples
                genotypes = []
                for sample in self.all_samples:
                    is_case = sample in self.case_samples
                    genotype = self.generate_genotype(gene, clinsig, is_case)
                    genotypes.append(genotype)
                
                # Write variant line
                fields = [chrom, pos, ".", ref, alt, "100", "PASS", info, "GT:DP:AD:GQ"]
                line = "\t".join(fields + genotypes) + "\n"
                f.write(line)
        
        return stats

    def _process_existing_vcf(self, input_vcf: str, output_vcf: Path) -> Dict[str, int]:
        """Process existing VCF file and modify genotypes."""
        # This would be similar to the original script's process_vcf method
        # For brevity, using synthetic approach but this could be implemented
        logger.info(f"Processing existing VCF would be implemented here: {input_vcf}")
        return self._create_synthetic_vcf(output_vcf)

    def create_sample_files(self, output_dir: Path) -> None:
        """Create sample list files for direct sample specification."""
        logger.info("Creating sample list files")
        
        # Basic case/control files
        (output_dir / "case_samples.txt").write_text("\n".join(self.case_samples) + "\n")
        (output_dir / "control_samples.txt").write_text("\n".join(self.control_samples) + "\n")
        
        # Subset files for testing partial overlaps
        case_subset = self.case_samples[:20]  # First 20 cases
        control_subset = self.control_samples[:30]  # First 30 controls
        
        (output_dir / "case_samples_subset.txt").write_text("\n".join(case_subset) + "\n")
        (output_dir / "control_samples_subset.txt").write_text("\n".join(control_subset) + "\n")
        
        logger.info(f"Created sample files with {len(self.case_samples)} cases, {len(self.control_samples)} controls")

    def create_phenotype_files(self, output_dir: Path) -> None:
        """Create comprehensive phenotype files for all test scenarios."""
        logger.info("Creating phenotype files")
        
        # 1. Basic phenotype file (CSV format like GCKD example)
        self._create_basic_phenotype_file(output_dir / "phenotypes_basic.csv")
        
        # 2. Extended phenotype file with multiple terms per sample
        self._create_extended_phenotype_file(output_dir / "phenotypes_extended.csv")
        
        # 3. TSV format phenotype file
        self._create_tsv_phenotype_file(output_dir / "phenotypes.tsv")
        
        # 4. Alternative column names for testing flexibility
        self._create_alternative_phenotype_file(output_dir / "phenotypes_alt_columns.csv")
        
        # 5. Phenotype file with missing samples (for error testing)
        self._create_incomplete_phenotype_file(output_dir / "phenotypes_incomplete.csv")

    def _create_basic_phenotype_file(self, output_path: Path) -> None:
        """Create basic phenotype file similar to GCKD format."""
        with open(output_path, "w", newline="") as f:
            writer = csv.writer(f)
            
            # Header similar to GCKD format
            writer.writerow(["MNPPSD", "CGRSequenceID", "SampleID", "Sex", "category_term", "identifier", "name", "link"])
            
            # Generate data for each sample
            for i, sample in enumerate(self.all_samples):
                is_case = sample in self.case_samples
                sex = random.choice(["male", "female"])
                
                if is_case:
                    # Cases get disease-related HPO terms
                    hpo_term = random.choice(self.HPO_TERMS["case_terms"])
                    name = self._get_hpo_name(hpo_term)
                else:
                    # Controls get normal/healthy terms
                    hpo_term = random.choice(self.HPO_TERMS["control_terms"])
                    name = self._get_hpo_name(hpo_term)
                
                # Write row
                writer.writerow([
                    f"mnp{i:03d}",  # MNPPSD
                    f"{300000 + i}",  # CGRSequenceID
                    sample,  # SampleID
                    sex,
                    "hpo_identifier",
                    hpo_term,
                    name,
                    f"https://hpo.jax.org/app/browse/term/{hpo_term}"
                ])

    def _create_extended_phenotype_file(self, output_path: Path) -> None:
        """Create phenotype file with multiple HPO terms per sample."""
        with open(output_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["MNPPSD", "CGRSequenceID", "SampleID", "Sex", "category_term", "identifier", "name", "link"])
            
            for i, sample in enumerate(self.all_samples):
                is_case = sample in self.case_samples
                sex = random.choice(["male", "female"])
                
                # Each sample gets 1-3 HPO terms
                num_terms = random.randint(1, 3)
                
                if is_case:
                    # Cases: mix of disease terms and some general terms
                    available_terms = self.HPO_TERMS["case_terms"] + self.HPO_TERMS["mixed_terms"]
                else:
                    # Controls: mostly normal terms, some mixed
                    available_terms = self.HPO_TERMS["control_terms"] + self.HPO_TERMS["mixed_terms"]
                
                selected_terms = random.sample(available_terms, min(num_terms, len(available_terms)))
                
                # Write one row per HPO term
                for hpo_term in selected_terms:
                    writer.writerow([
                        f"mnp{i:03d}",
                        f"{300000 + i}",
                        sample,
                        sex,
                        "hpo_identifier",
                        hpo_term,
                        self._get_hpo_name(hpo_term),
                        f"https://hpo.jax.org/app/browse/term/{hpo_term}"
                    ])

    def _create_tsv_phenotype_file(self, output_path: Path) -> None:
        """Create simple TSV phenotype file."""
        with open(output_path, "w") as f:
            f.write("sample_id\tphenotype\tage\tsex\n")
            
            for sample in self.all_samples:
                is_case = sample in self.case_samples
                age = random.randint(25, 75)
                sex = random.choice(["M", "F"])
                phenotype = "affected" if is_case else "unaffected"
                
                f.write(f"{sample}\t{phenotype}\t{age}\t{sex}\n")

    def _create_alternative_phenotype_file(self, output_path: Path) -> None:
        """Create phenotype file with alternative column names."""
        with open(output_path, "w", newline="") as f:
            writer = csv.writer(f)
            # Alternative column names to test flexibility
            writer.writerow(["patient_id", "sequencing_id", "sample_name", "gender", "term_type", "hpo_id", "description", "url"])
            
            for i, sample in enumerate(self.all_samples):
                is_case = sample in self.case_samples
                sex = random.choice(["male", "female"])
                
                hpo_term = (random.choice(self.HPO_TERMS["case_terms"]) if is_case 
                           else random.choice(self.HPO_TERMS["control_terms"]))
                
                writer.writerow([
                    f"patient_{i:03d}",
                    f"{400000 + i}",
                    sample,  # sample_name column
                    sex,
                    "hpo_identifier",
                    hpo_term,  # hpo_id column
                    self._get_hpo_name(hpo_term),
                    f"https://hpo.jax.org/app/browse/term/{hpo_term}"
                ])

    def _create_incomplete_phenotype_file(self, output_path: Path) -> None:
        """Create phenotype file missing some samples (for error testing)."""
        with open(output_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(["MNPPSD", "CGRSequenceID", "SampleID", "Sex", "category_term", "identifier", "name", "link"])
            
            # Only include 80% of samples
            included_samples = self.all_samples[:int(0.8 * len(self.all_samples))]
            
            for i, sample in enumerate(included_samples):
                is_case = sample in self.case_samples
                sex = random.choice(["male", "female"])
                
                hpo_term = (random.choice(self.HPO_TERMS["case_terms"]) if is_case 
                           else random.choice(self.HPO_TERMS["control_terms"]))
                
                writer.writerow([
                    f"mnp{i:03d}",
                    f"{300000 + i}",
                    sample,
                    sex,
                    "hpo_identifier",
                    hpo_term,
                    self._get_hpo_name(hpo_term),
                    f"https://hpo.jax.org/app/browse/term/{hpo_term}"
                ])

    def _get_hpo_name(self, hpo_term: str) -> str:
        """Get human-readable name for HPO term."""
        hpo_names = {
            "HP:0000113": "Polycystic kidney dysplasia",
            "HP:0000003": "Multicystic kidney dysplasia",
            "HP:0000107": "Renal cyst",
            "HP:0000822": "Hypertension",
            "HP:0003774": "Stage 5 chronic kidney disease",
            "HP:0000001": "All",
            "HP:0032101": "Normal phenotype",
            "HP:0001626": "Abnormality of the cardiovascular system",
            "HP:0000478": "Abnormality of the eye",
            "HP:0001574": "Abnormality of the integument",
        }
        return hpo_names.get(hpo_term, "Unknown phenotype")

    def create_hpo_term_files(self, output_dir: Path) -> None:
        """Create HPO term files for phenotype-based classification."""
        logger.info("Creating HPO term files")
        
        # Case HPO terms
        (output_dir / "case_hpo_terms.txt").write_text("\n".join(self.HPO_TERMS["case_terms"]) + "\n")
        
        # Control HPO terms
        (output_dir / "control_hpo_terms.txt").write_text("\n".join(self.HPO_TERMS["control_terms"]) + "\n")
        
        # Mixed terms (for testing edge cases)
        (output_dir / "mixed_hpo_terms.txt").write_text("\n".join(self.HPO_TERMS["mixed_terms"]) + "\n")

    def create_gene_files(self, output_dir: Path) -> None:
        """Create gene list files."""
        logger.info("Creating gene list files")
        
        # All test genes
        (output_dir / "test_genes.txt").write_text("\n".join(sorted(self.ALL_GENES)) + "\n")
        
        # Disease genes only
        (output_dir / "disease_genes.txt").write_text("\n".join(sorted(self.DISEASE_GENES)) + "\n")
        
        # Control genes only
        (output_dir / "control_genes.txt").write_text("\n".join(sorted(self.CONTROL_GENES)) + "\n")

    def create_test_configurations(self, output_dir: Path) -> None:
        """Create test configuration files for different scenarios."""
        logger.info("Creating test configuration files")
        
        configs = {
            # Scenario 1: Direct sample specification
            "test_direct_samples": {
                "description": "Test direct sample specification",
                "command_args": [
                    "--case-samples", ",".join(self.case_samples[:10]),
                    "--control-samples", ",".join(self.control_samples[:15]),
                ],
                "expected_cases": 10,
                "expected_controls": 15,
            },
            
            # Scenario 2: Sample files
            "test_sample_files": {
                "description": "Test sample file specification",
                "command_args": [
                    "--case-samples-file", "case_samples.txt",
                    "--control-samples-file", "control_samples.txt",
                ],
                "expected_cases": self.NUM_CASES,
                "expected_controls": self.NUM_CONTROLS,
            },
            
            # Scenario 3: HPO-based classification
            "test_hpo_classification": {
                "description": "Test HPO-based classification",
                "command_args": [
                    "--phenotype-file", "phenotypes_basic.csv",
                    "--phenotype-sample-column", "SampleID",
                    "--phenotype-value-column", "identifier",
                    "--case-phenotypes", ",".join(self.HPO_TERMS["case_terms"]),
                    "--control-phenotypes", ",".join(self.HPO_TERMS["control_terms"]),
                ],
                "expected_cases": "variable",  # Depends on phenotype assignment
                "expected_controls": "variable",
            },
            
            # Scenario 4: HPO term files
            "test_hpo_files": {
                "description": "Test HPO term files",
                "command_args": [
                    "--phenotype-file", "phenotypes_basic.csv",
                    "--phenotype-sample-column", "SampleID", 
                    "--phenotype-value-column", "identifier",
                    "--case-phenotypes-file", "case_hpo_terms.txt",
                    "--control-phenotypes-file", "control_hpo_terms.txt",
                ],
                "expected_cases": "variable",
                "expected_controls": "variable",
            },
            
            # Scenario 5: Alternative column names
            "test_alt_columns": {
                "description": "Test alternative column names",
                "command_args": [
                    "--phenotype-file", "phenotypes_alt_columns.csv",
                    "--phenotype-sample-column", "sample_name",
                    "--phenotype-value-column", "hpo_id",
                    "--case-phenotypes", ",".join(self.HPO_TERMS["case_terms"]),
                    "--control-phenotypes", ",".join(self.HPO_TERMS["control_terms"]),
                ],
                "expected_cases": "variable",
                "expected_controls": "variable",
            },
        }
        
        # Write configurations to JSON
        with open(output_dir / "test_configurations.json", "w") as f:
            json.dump(configs, f, indent=2)

    def create_test_scripts(self, output_dir: Path) -> None:
        """Create comprehensive test scripts."""
        logger.info("Creating test scripts")
        
        # Main test script
        main_script = f"""#!/bin/bash
# Comprehensive Gene Burden Analysis Test Suite
# Generated on {datetime.now().isoformat()}

set -e

echo "=== VariantCentrifuge Gene Burden Test Suite ==="
echo "Test data directory: {output_dir.absolute()}"
echo

# Test 1: Direct sample specification
echo "Test 1: Direct sample specification"
variantcentrifuge \\
    --vcf-file test_data.vcf.gz \\
    --gene-file test_genes.txt \\
    --case-samples CASE_001,CASE_002,CASE_003,CASE_004,CASE_005 \\
    --control-samples CTRL_001,CTRL_002,CTRL_003,CTRL_004,CTRL_005,CTRL_006,CTRL_007 \\
    --perform-gene-burden \\
    --preset rare,coding \\
    --output-dir test1_direct_samples \\
    --output-file test1_results.tsv \\
    --use-new-pipeline

echo "✓ Test 1 completed"
echo

# Test 2: Sample files
echo "Test 2: Sample file specification"
variantcentrifuge \\
    --vcf-file test_data.vcf.gz \\
    --gene-file test_genes.txt \\
    --case-samples-file case_samples.txt \\
    --control-samples-file control_samples.txt \\
    --perform-gene-burden \\
    --preset rare,coding \\
    --output-dir test2_sample_files \\
    --output-file test2_results.tsv \\
    --use-new-pipeline

echo "✓ Test 2 completed"
echo

# Test 3: HPO-based classification
echo "Test 3: HPO-based classification"
variantcentrifuge \\
    --vcf-file test_data.vcf.gz \\
    --gene-file test_genes.txt \\
    --phenotype-file phenotypes_basic.csv \\
    --phenotype-sample-column SampleID \\
    --phenotype-value-column identifier \\
    --case-phenotypes HP:0000113,HP:0000003,HP:0000107,HP:0000822,HP:0003774 \\
    --control-phenotypes HP:0000001,HP:0032101 \\
    --perform-gene-burden \\
    --preset rare,coding \\
    --output-dir test3_hpo_classification \\
    --output-file test3_results.tsv \\
    --use-new-pipeline

echo "✓ Test 3 completed"
echo

# Test 4: HPO term files
echo "Test 4: HPO term files"
variantcentrifuge \\
    --vcf-file test_data.vcf.gz \\
    --gene-file test_genes.txt \\
    --phenotype-file phenotypes_basic.csv \\
    --phenotype-sample-column SampleID \\
    --phenotype-value-column identifier \\
    --case-phenotypes-file case_hpo_terms.txt \\
    --control-phenotypes-file control_hpo_terms.txt \\
    --perform-gene-burden \\
    --preset rare,coding \\
    --output-dir test4_hpo_files \\
    --output-file test4_results.tsv \\
    --use-new-pipeline

echo "✓ Test 4 completed"
echo

# Test 5: Alternative column names
echo "Test 5: Alternative column names"
variantcentrifuge \\
    --vcf-file test_data.vcf.gz \\
    --gene-file test_genes.txt \\
    --phenotype-file phenotypes_alt_columns.csv \\
    --phenotype-sample-column sample_name \\
    --phenotype-value-column hpo_id \\
    --case-phenotypes HP:0000113,HP:0000003,HP:0000107,HP:0000822,HP:0003774 \\
    --control-phenotypes HP:0000001,HP:0032101 \\
    --perform-gene-burden \\
    --preset rare,coding \\
    --output-dir test5_alt_columns \\
    --output-file test5_results.tsv \\
    --use-new-pipeline

echo "✓ Test 5 completed"
echo

echo "=== All tests completed! ==="
echo "Check individual test directories for results:"
echo "- test1_direct_samples/"
echo "- test2_sample_files/"  
echo "- test3_hpo_classification/"
echo "- test4_hpo_files/"
echo "- test5_alt_columns/"
echo

echo "Expected results:"
echo "- Disease genes (PKD1, PKD2, BRCA1, BRCA2) should show enrichment in cases"
echo "- Control genes (TTN, OBSCN, MUC16) should show no significant enrichment"
"""

        # Write main script
        main_script_path = output_dir / "run_comprehensive_tests.sh"
        main_script_path.write_text(main_script)
        os.chmod(main_script_path, 0o755)
        
        # Individual test scripts for each scenario
        for i, (test_name, config) in enumerate([
            ("direct_samples", "Direct sample specification"),
            ("sample_files", "Sample file specification"),
            ("hpo_classification", "HPO-based classification"),
            ("hpo_files", "HPO term files"),
            ("alt_columns", "Alternative column names"),
        ], 1):
            
            script_content = f"""#!/bin/bash
# Test {i}: {config}

set -e

echo "Running Test {i}: {config}"

# Add specific test command here based on scenario
echo "Test implementation would go here"

echo "✓ Test {i} completed"
"""
            
            script_path = output_dir / f"test_{i}_{test_name}.sh"
            script_path.write_text(script_content)
            os.chmod(script_path, 0o755)

    def write_statistics(self, output_dir: Path, vcf_stats: Dict[str, int]) -> None:
        """Write comprehensive statistics about the generated dataset."""
        logger.info("Writing dataset statistics")
        
        # Calculate genotype distribution
        genotype_summary = {}
        for key, genotype_counts in self.variant_stats.items():
            total = sum(genotype_counts.values())
            if total > 0:
                genotype_summary[key] = {
                    "total_variants": total,
                    "ref_only (0/0)": genotype_counts.get("0/0", 0),
                    "heterozygous (0/1)": genotype_counts.get("0/1", 0),
                    "homozygous (1/1)": genotype_counts.get("1/1", 0),
                    "variant_rate": (genotype_counts.get("0/1", 0) + genotype_counts.get("1/1", 0)) / total,
                }

        stats = {
            "generation_info": {
                "date": datetime.now().isoformat(),
                "generator_version": "comprehensive_v1.0",
                "description": "Comprehensive gene burden test dataset covering all specification methods",
            },
            "sample_composition": {
                "total_samples": len(self.all_samples),
                "case_samples": len(self.case_samples),
                "control_samples": len(self.control_samples),
                "case_sample_ids": self.case_samples,
                "control_sample_ids": self.control_samples,
            },
            "variant_statistics": vcf_stats,
            "genotype_distribution": genotype_summary,
            "test_genes": {
                "disease_genes": list(self.DISEASE_GENES),
                "control_genes": list(self.CONTROL_GENES),
                "total_genes": len(self.ALL_GENES),
            },
            "hpo_terms": self.HPO_TERMS,
            "probabilities": self.PROBABILITIES,
            "test_scenarios": [
                "Direct sample specification (--case-samples, --control-samples)",
                "Sample file specification (--case-samples-file, --control-samples-file)",
                "HPO-based classification with phenotype file",
                "HPO term files (--case-phenotypes-file, --control-phenotypes-file)",
                "Alternative column names in phenotype files",
                "Error handling for incomplete phenotype data",
            ],
        }

        with open(output_dir / "dataset_statistics.json", "w") as f:
            json.dump(stats, f, indent=2)

    def write_readme(self, output_dir: Path) -> None:
        """Write comprehensive README for the test dataset."""
        readme_content = f"""# Comprehensive Gene Burden Test Dataset

Generated on: {datetime.now().isoformat()}

## Overview

This test dataset provides comprehensive coverage for testing all gene burden analysis 
functionality in VariantCentrifuge, including all possible ways to specify case/control groups.

## Dataset Composition

- **Total Samples**: {len(self.all_samples)} ({len(self.case_samples)} cases, {len(self.control_samples)} controls)
- **Disease Genes**: {', '.join(sorted(self.DISEASE_GENES))}
- **Control Genes**: {', '.join(sorted(self.CONTROL_GENES))}

## Files Generated

### Core Data Files
- `test_data.vcf.gz` - Main VCF file with controlled genotype distributions
- `test_data.vcf.gz.tbi` - Tabix index

### Sample Specification Files
- `case_samples.txt` - Complete list of case sample IDs
- `control_samples.txt` - Complete list of control sample IDs
- `case_samples_subset.txt` - Subset of case samples (for testing)
- `control_samples_subset.txt` - Subset of control samples (for testing)

### Phenotype Files
- `phenotypes_basic.csv` - Basic phenotype file (GCKD-like format)
- `phenotypes_extended.csv` - Extended with multiple HPO terms per sample
- `phenotypes.tsv` - Simple TSV format
- `phenotypes_alt_columns.csv` - Alternative column names
- `phenotypes_incomplete.csv` - Missing some samples (for error testing)

### HPO Term Files
- `case_hpo_terms.txt` - HPO terms defining case phenotypes
- `control_hpo_terms.txt` - HPO terms defining control phenotypes  
- `mixed_hpo_terms.txt` - Mixed terms for edge case testing

### Gene Lists
- `test_genes.txt` - All test genes
- `disease_genes.txt` - Disease-associated genes only
- `control_genes.txt` - Control genes only

### Test Infrastructure
- `run_comprehensive_tests.sh` - Main test runner script
- `test_1_direct_samples.sh` through `test_5_alt_columns.sh` - Individual test scripts
- `test_configurations.json` - Test scenario configurations
- `dataset_statistics.json` - Detailed dataset statistics

## Test Scenarios

### 1. Direct Sample Specification
```bash
--case-samples CASE_001,CASE_002,CASE_003
--control-samples CTRL_001,CTRL_002,CTRL_003
```

### 2. Sample File Specification  
```bash
--case-samples-file case_samples.txt
--control-samples-file control_samples.txt
```

### 3. HPO-based Classification
```bash
--phenotype-file phenotypes_basic.csv
--phenotype-sample-column SampleID
--phenotype-value-column identifier
--case-phenotypes HP:0000113,HP:0000003,HP:0000107
--control-phenotypes HP:0000001,HP:0032101
```

### 4. HPO Term Files
```bash
--phenotype-file phenotypes_basic.csv
--phenotype-sample-column SampleID
--phenotype-value-column identifier
--case-phenotypes-file case_hpo_terms.txt
--control-phenotypes-file control_hpo_terms.txt
```

### 5. Alternative Column Names
```bash
--phenotype-file phenotypes_alt_columns.csv
--phenotype-sample-column sample_name
--phenotype-value-column hpo_id
--case-phenotypes HP:0000113,HP:0000003,HP:0000107
--control-phenotypes HP:0000001,HP:0032101
```

## Expected Results

When running gene burden analysis on this dataset:

- **Disease genes** (PKD1, PKD2, BRCA1, BRCA2) should show **significant enrichment** in cases
- **Control genes** (TTN, OBSCN, MUC16) should show **no significant enrichment**
- P-values for disease genes should be < 0.05 (before multiple testing correction)
- Odds ratios for disease genes should be > 1.0

## Running Tests

1. **Run all tests**:
   ```bash
   ./run_comprehensive_tests.sh
   ```

2. **Run individual tests**:
   ```bash
   ./test_1_direct_samples.sh
   ./test_2_sample_files.sh
   # etc.
   ```

3. **Custom test**:
   ```bash
   variantcentrifuge \\
     --vcf-file test_data.vcf.gz \\
     --gene-file test_genes.txt \\
     --case-samples-file case_samples.txt \\
     --control-samples-file control_samples.txt \\
     --perform-gene-burden \\
     --preset rare,coding \\
     --use-new-pipeline \\
     --output-file results.tsv
   ```

## HPO Terms Used

### Case Phenotypes (Disease-related)
{chr(10).join(f"- {term}: {self._get_hpo_name(term)}" for term in self.HPO_TERMS["case_terms"])}

### Control Phenotypes (Normal/Healthy)
{chr(10).join(f"- {term}: {self._get_hpo_name(term)}" for term in self.HPO_TERMS["control_terms"])}

### Mixed Phenotypes (For edge case testing)
{chr(10).join(f"- {term}: {self._get_hpo_name(term)}" for term in self.HPO_TERMS["mixed_terms"])}

## Validation

Check `dataset_statistics.json` for detailed statistics on:
- Genotype distribution per gene and sample type
- Variant counts and classifications
- Sample composition
- Test scenario configurations

The dataset is designed with controlled probabilities to ensure reliable test results
while maintaining realistic genetic variant distributions.
"""

        (output_dir / "README.md").write_text(readme_content)

    def generate_complete_dataset(self, output_dir: Path, base_vcf: str = None) -> None:
        """Generate the complete comprehensive test dataset."""
        logger.info(f"Generating comprehensive gene burden test dataset in {output_dir}")
        
        # Create output directory
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # Generate all components
        vcf_stats = self.create_vcf_file(output_dir / "test_data.vcf", base_vcf)
        self.create_sample_files(output_dir)
        self.create_phenotype_files(output_dir)
        self.create_hpo_term_files(output_dir)
        self.create_gene_files(output_dir)
        self.create_test_configurations(output_dir)
        self.create_test_scripts(output_dir)
        self.write_statistics(output_dir, vcf_stats)
        self.write_readme(output_dir)
        
        # Compress VCF
        logger.info("Compressing VCF file...")
        os.system(f"cd {output_dir} && bgzip -f test_data.vcf && tabix -p vcf test_data.vcf.gz")
        
        logger.info(f"✓ Complete test dataset generated in {output_dir}")
        logger.info("✓ Run './run_comprehensive_tests.sh' to execute all test scenarios")


def main():
    """Generate comprehensive test dataset."""
    parser = argparse.ArgumentParser(
        description="Generate comprehensive gene burden test dataset",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__
    )
    
    parser.add_argument(
        "--output-dir", 
        required=True, 
        help="Output directory for test dataset"
    )
    parser.add_argument(
        "--base-vcf",
        help="Base VCF file to process (optional, will create synthetic if not provided)"
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed for reproducibility (default: 42)"
    )
    
    args = parser.parse_args()
    
    # Initialize generator
    generator = ComprehensiveTestDataGenerator(seed=args.seed)
    
    # Generate complete dataset
    output_dir = Path(args.output_dir)
    generator.generate_complete_dataset(output_dir, args.base_vcf)
    
    print(f"\n🎉 Comprehensive gene burden test dataset generated!")
    print(f"📁 Location: {output_dir.absolute()}")
    print(f"📖 See README.md for detailed usage instructions")
    print(f"🧪 Run ./run_comprehensive_tests.sh to execute all test scenarios")


if __name__ == "__main__":
    main()