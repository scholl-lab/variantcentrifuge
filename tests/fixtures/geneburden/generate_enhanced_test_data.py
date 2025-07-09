#!/usr/bin/env python3
"""
Enhanced test data generator with realistic annotation sampling for gene burden analysis.

This script creates realistic test datasets by:
1. Sampling authentic SnpEff annotations from anonymized local templates
2. Mapping limited source genes to diverse target genes
3. Preserving complex annotation structures while ensuring controlled test outcomes
4. Maintaining comprehensive test coverage for all gene burden scenarios

Key Features:
- Real annotation pattern sampling with authentic dbNSFP scores
- Gene mapping with effect-aware annotation transfer
- Modular annotation sampling classes for reusability
- Controlled genotype distributions for reliable test outcomes
- Comprehensive test scenario coverage
- Self-contained with local anonymized annotation templates

Author: Enhanced by Claude Code
"""

import argparse
import csv
import json
import logging
import os
import random
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional, Tuple

# Configure logging
logging.basicConfig(level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")
logger = logging.getLogger(__name__)


@dataclass
class AnnotationTemplate:
    """Template for variant annotations extracted from real data."""

    effect: str  # e.g., "missense_variant"
    impact: str  # HIGH, MODERATE, LOW, MODIFIER
    gene_biotype: str  # protein_coding, pseudogene, etc.
    transcript_biotype: str  # protein_coding, etc.
    ann_fields: str  # Full ANN field content
    dbNSFP_fields: Dict[str, str]  # All dbNSFP prediction scores
    other_fields: Dict[str, str]  # Other INFO fields
    consequence_rank: int  # Ranking of consequence severity


class AnnotationSampler:
    """Samples and categorizes real annotation patterns from test data."""

    def __init__(self, source_file: str = None):
        """Initialize with source annotation file."""
        # Default to local anonymized annotation templates
        if source_file is None:
            script_dir = Path(__file__).parent
            source_file = str(script_dir / "test_data" / "annotation_templates.txt")

        self.source_file = source_file
        self.templates = defaultdict(list)  # effect_type -> [AnnotationTemplate]
        self.gene_annotations = defaultdict(list)  # source_gene -> [AnnotationTemplate]
        self.consequence_ranking = {
            # High impact
            "stop_gained": 10,
            "frameshift_variant": 10,
            "stop_lost": 10,
            "start_lost": 10,
            "splice_acceptor_variant": 10,
            "splice_donor_variant": 10,
            # Moderate impact
            "missense_variant": 7,
            "inframe_deletion": 6,
            "inframe_insertion": 6,
            "splice_region_variant": 5,
            # Low impact
            "synonymous_variant": 3,
            "stop_retained_variant": 3,
            # Modifier impact
            "intron_variant": 1,
            "upstream_gene_variant": 1,
            "downstream_gene_variant": 1,
            "intergenic_variant": 0,
        }

        logger.info(f"Initializing annotation sampler with {source_file}")
        self._parse_annotations()

    def _parse_annotations(self) -> None:
        """Parse annotations from source file and categorize them."""
        logger.info("Parsing real annotations from source file...")

        parsed_count = 0

        try:
            with open(self.source_file, "r") as f:
                for line_num, line in enumerate(f, 1):
                    if line.startswith("#") or not line.strip():
                        continue

                    try:
                        # Handle both VCF format (8+ fields) and INFO-only format
                        line = line.strip()
                        if "\t" in line:
                            fields = line.split("\t")
                            if len(fields) >= 8:
                                info = fields[7]  # Standard VCF INFO field
                            else:
                                continue
                        else:
                            # INFO field only format
                            info = line

                        # Extract annotation template
                        template = self._extract_annotation_template(info)
                        if template:
                            # Categorize by effect type
                            self.templates[template.effect].append(template)

                            # Extract gene from ANN field
                            gene = self._extract_gene_from_ann(template.ann_fields)
                            if gene:
                                self.gene_annotations[gene].append(template)

                            parsed_count += 1

                    except Exception as e:
                        logger.warning(f"Error parsing line {line_num}: {e}")
                        continue

        except FileNotFoundError:
            logger.error(f"Source file not found: {self.source_file}")
            raise

        logger.info(f"Parsed {parsed_count} annotation templates")
        logger.info(f"Found {len(self.templates)} unique effect types")
        logger.info(f"Found {len(self.gene_annotations)} source genes")

        # Log effect distribution
        for effect, templates in self.templates.items():
            logger.debug(f"  {effect}: {len(templates)} templates")

    def _extract_annotation_template(self, info_field: str) -> Optional[AnnotationTemplate]:
        """Extract annotation template from INFO field."""
        try:
            # Extract ANN field (handle both ANN= and mangled MMLEMLEANN= formats)
            ann_match = re.search(r"(?:ANN|MMLEMLEANN|[A-Z]*ANN)=([^;]+)", info_field)
            if not ann_match:
                return None

            ann_value = ann_match.group(1)

            # Parse first annotation (primary transcript)
            ann_parts = ann_value.split("|")
            if len(ann_parts) < 10:
                return None

            alt, effect, impact, gene, gene_id, feature_type, feature_id, transcript_biotype = (
                ann_parts[:8]
            )

            # Extract ALL annotation fields (dbNSFP, splice, ClinVar, HGMD, etc.)
            dbNSFP_fields = {}
            other_fields = {}

            # Extract dbNSFP fields
            for match in re.finditer(r"(dbNSFP_[^=]+)=([^;]+)", info_field):
                key, value = match.groups()
                dbNSFP_fields[key] = value

            # Extract splice prediction fields
            for match in re.finditer(r"(splice_[^=]+)=([^;]+)", info_field):
                key, value = match.groups()
                other_fields[key] = value

            # Extract ClinVar and HGMD fields
            for pattern in [
                r"(ClinVar_[^=]+)=([^;]+)",
                r"(hgmd_[^=]+)=([^;]+)",
                r"(CADD_phred)=([^;]+)",
                r"(gnomAD_[^=]+)=([^;]+)",
                r"(REVEL_[^=]+)=([^;]+)",
                r"(SIFT[^=]+)=([^;]+)",
                r"(Polyphen2_[^=]+)=([^;]+)",
                r"(FATHMM_[^=]+)=([^;]+)",
            ]:
                for match in re.finditer(pattern, info_field):
                    key, value = match.groups()
                    other_fields[key] = value

            # Rank consequence severity
            consequence_rank = self.consequence_ranking.get(effect, 0)

            template = AnnotationTemplate(
                effect=effect,
                impact=impact,
                gene_biotype=feature_type,
                transcript_biotype=transcript_biotype,
                ann_fields=ann_value,
                dbNSFP_fields=dbNSFP_fields,
                other_fields=other_fields,
                consequence_rank=consequence_rank,
            )

            return template

        except Exception as e:
            logger.debug(f"Error extracting annotation template: {e}")
            return None

    def _extract_gene_from_ann(self, ann_field: str) -> Optional[str]:
        """Extract gene name from ANN field."""
        try:
            parts = ann_field.split("|")
            if len(parts) > 3:
                return parts[3]  # Gene name is in position 3
        except Exception:
            pass
        return None

    def sample_annotation(
        self, target_gene: str, effect_preference: str = None, impact_preference: str = None
    ) -> Optional[AnnotationTemplate]:
        """Sample an appropriate annotation template."""

        # Get candidates based on preferences
        candidates = []

        if effect_preference and effect_preference in self.templates:
            candidates = self.templates[effect_preference]
        elif impact_preference:
            # Find templates with matching impact
            for effect, templates in self.templates.items():
                candidates.extend([t for t in templates if t.impact == impact_preference])
        else:
            # Use all available templates
            for templates in self.templates.values():
                candidates.extend(templates)

        if not candidates:
            logger.warning(f"No annotation candidates found for {target_gene}")
            return None

        # Select template and adapt for target gene
        template = random.choice(candidates)

        # Create adapted template
        adapted_ann = self._adapt_annotation_for_gene(template.ann_fields, target_gene)

        adapted_template = AnnotationTemplate(
            effect=template.effect,
            impact=template.impact,
            gene_biotype=template.gene_biotype,
            transcript_biotype=template.transcript_biotype,
            ann_fields=adapted_ann,
            dbNSFP_fields=template.dbNSFP_fields.copy(),
            other_fields=template.other_fields.copy(),
            consequence_rank=template.consequence_rank,
        )

        return adapted_template

    def _adapt_annotation_for_gene(self, ann_field: str, target_gene: str) -> str:
        """Adapt annotation field for target gene."""
        parts = ann_field.split("|")
        if len(parts) > 3:
            parts[3] = target_gene  # Replace gene name
            parts[4] = target_gene  # Replace gene ID
        return "|".join(parts)

    def get_effect_types(self) -> List[str]:
        """Get all available effect types."""
        return list(self.templates.keys())

    def get_high_impact_effects(self) -> List[str]:
        """Get high impact effect types."""
        return [
            effect
            for effect, templates in self.templates.items()
            if templates and templates[0].impact == "HIGH"
        ]

    def get_moderate_impact_effects(self) -> List[str]:
        """Get moderate impact effect types."""
        return [
            effect
            for effect, templates in self.templates.items()
            if templates and templates[0].impact == "MODERATE"
        ]


class EnhancedTestDataGenerator:
    """Enhanced test data generator using real annotation sampling."""

    # Test genes with real hg19 coordinates
    DISEASE_GENES = {"PKD1", "PKD2", "BRCA1", "BRCA2"}
    CONTROL_GENES = {"TTN", "OBSCN", "MUC16"}
    ALL_GENES = DISEASE_GENES | CONTROL_GENES

    # Real hg19 gene coordinates (chrom, start, end, strand)
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
        "case_disease_benign": 0.25,  # 25% cases have benign disease gene variants
        "control_disease_pathogenic": 0.10,  # 10% controls have pathogenic disease variants
        "control_disease_benign": 0.15,  # 15% controls have benign disease variants
        "any_control_gene": 0.40,  # 40% all samples have control gene variants
    }

    HOMO_PROB = 0.20  # Probability of homozygous when variant is present

    def __init__(self, annotation_source: str = None, seed: int = 42):
        """Initialize generator with annotation sampler."""
        random.seed(seed)

        # Initialize annotation sampler (defaults to local anonymized templates)
        self.annotation_sampler = AnnotationSampler(annotation_source)

        # Generate sample IDs with realistic naming
        self.case_samples = [f"CASE_{i:03d}" for i in range(1, self.NUM_CASES + 1)]
        self.control_samples = [f"CTRL_{i:03d}" for i in range(1, self.NUM_CONTROLS + 1)]
        self.all_samples = self.case_samples + self.control_samples

        # Track statistics
        self.variant_stats = defaultdict(lambda: defaultdict(int))

        logger.info(
            f"Initialized enhanced generator with {len(self.case_samples)} cases, "
            f"{len(self.control_samples)} controls"
        )
        logger.info(f"Annotation templates loaded: {len(self.annotation_sampler.templates)}")

    def generate_realistic_variants_with_annotations(
        self,
    ) -> List[Tuple[str, str, str, str, str, str]]:
        """Generate variants with realistic annotations sampled from real data."""
        variants = []

        # Predefined alleles for realistic substitutions
        alleles = ["A", "T", "G", "C"]

        for gene in self.ALL_GENES:
            if gene not in self.GENE_COORDINATES:
                logger.warning(f"No coordinates found for gene {gene}, skipping")
                continue

            chrom, start, end, strand = self.GENE_COORDINATES[gene]
            gene_length = end - start

            # Determine variant characteristics based on gene type
            if gene in self.DISEASE_GENES:
                # Disease genes: more variants with diverse impacts
                num_variants = 4 if gene_length < 100000 else 5

                # Create mix of high/moderate impact variants for disease genes
                variant_effects = []

                # Add high impact variants (pathogenic)
                high_impact_effects = self.annotation_sampler.get_high_impact_effects()
                if high_impact_effects:
                    variant_effects.extend(random.choices(high_impact_effects, k=2))

                # Add moderate impact variants
                moderate_impact_effects = self.annotation_sampler.get_moderate_impact_effects()
                if moderate_impact_effects:
                    variant_effects.extend(
                        random.choices(moderate_impact_effects, k=num_variants - 2)
                    )

                # Fill remaining with any available effects
                while len(variant_effects) < num_variants:
                    variant_effects.append(
                        random.choice(self.annotation_sampler.get_effect_types())
                    )

            else:
                # Control genes: mostly benign/moderate effects
                num_variants = 3
                variant_effects = random.choices(
                    self.annotation_sampler.get_moderate_impact_effects()
                    or self.annotation_sampler.get_effect_types(),
                    k=num_variants,
                )

            # Generate variants within gene boundaries
            for i in range(num_variants):
                # Sample position within gene (avoid very start/end for realism)
                margin = min(gene_length // 10, 5000)
                pos_start = start + margin
                pos_end = end - margin
                position = random.randint(pos_start, pos_end)

                # Generate realistic alleles
                ref_allele = random.choice(alleles)
                alt_allele = random.choice([a for a in alleles if a != ref_allele])

                # Get realistic annotation for this effect and gene
                effect_type = (
                    variant_effects[i]
                    if i < len(variant_effects)
                    else random.choice(self.annotation_sampler.get_effect_types())
                )
                annotation = self.annotation_sampler.sample_annotation(
                    target_gene=gene, effect_preference=effect_type
                )

                if annotation:
                    variants.append(
                        (chrom, str(position), ref_allele, alt_allele, gene, annotation)
                    )
                else:
                    logger.warning(
                        f"Could not generate annotation for {gene} with effect {effect_type}"
                    )

        # Sort variants by chromosome and position
        def sort_key(variant):
            chrom, pos, ref, alt, gene, annotation = variant
            try:
                chrom_num = (
                    int(chrom)
                    if chrom.isdigit()
                    else (23 if chrom == "X" else 24 if chrom == "Y" else 25)
                )
                return (chrom_num, int(pos))
            except Exception:
                return (999, 0)

        variants.sort(key=sort_key)
        logger.info(
            f"Generated {len(variants)} variants with real annotations across {len(self.ALL_GENES)} genes"
        )

        return variants

    def create_info_field(
        self, annotation: AnnotationTemplate, gene: str, genotypes: List[str]
    ) -> str:
        """Create realistic INFO field from annotation template."""
        info_parts = []

        # Calculate AC (allele count) from genotypes
        ac = 0
        for genotype in genotypes:
            gt = genotype.split(":")[0]  # Extract GT field
            if gt == "0/1":
                ac += 1
            elif gt == "1/1":
                ac += 2

        # Add core fields
        info_parts.append(f"AC={ac}")

        # Properly escape ANN field content to avoid parsing issues
        ann_content = annotation.ann_fields.replace(";", "%3B").replace("=", "%3D")
        info_parts.append(f"ANN={ann_content}")
        info_parts.extend(["SNP", "VARTYPE=SNP"])

        # Add dbNSFP fields (skip any that might cause issues)
        for key, value in annotation.dbNSFP_fields.items():
            # Skip keys that contain problematic characters or are too long
            if "|" in key or "&" in key or len(key) > 50:
                continue
            clean_key = key.replace("|", "_").replace("&", "_").replace("=", "_")
            clean_value = str(value).replace(";", ",").replace("=", ":")
            info_parts.append(f"{clean_key}={clean_value}")

        # Add other annotation fields (also filtered)
        for key, value in annotation.other_fields.items():
            # Skip keys that contain problematic characters or are too long
            if "|" in key or "&" in key or len(key) > 50:
                continue
            clean_key = key.replace("|", "_").replace("&", "_").replace("=", "_")
            clean_value = str(value).replace(";", ",").replace("=", ":")
            info_parts.append(f"{clean_key}={clean_value}")

        # Add default values for commonly expected fields if missing
        expected_fields = {
            "splice_dbscSNV_rf_score": "0.0",
            "splice_dbscSNV_ada_score": "0.0",
            "splice_spidex_dpsi_zscore": "0.0",
            "dbNSFP_gnomAD_exomes_AC": "0",
            "dbNSFP_gnomAD_genomes_AC": "0",
            "dbNSFP_ALFA_Total_AC": "0",
            "dbNSFP_clinvar_clnsig": "uncertain",
            "ClinVar_CLNSIG": "uncertain",
            "hgmd_CLASS": "unknown",
        }

        # Get all existing field names
        existing_fields = set()
        for key in annotation.dbNSFP_fields.keys():
            if not ("|" in key or "&" in key or len(key) > 50):
                existing_fields.add(key.replace("|", "_").replace("&", "_").replace("=", "_"))
        for key in annotation.other_fields.keys():
            if not ("|" in key or "&" in key or len(key) > 50):
                existing_fields.add(key.replace("|", "_").replace("&", "_").replace("=", "_"))

        # Add missing expected fields with defaults
        for field, default_value in expected_fields.items():
            if field not in existing_fields:
                info_parts.append(f"{field}={default_value}")

        return ";".join(info_parts)

    def is_pathogenic_annotation(self, annotation: AnnotationTemplate) -> bool:
        """Determine if annotation suggests pathogenic variant."""
        # High impact variants are considered pathogenic
        if annotation.impact == "HIGH":
            return True

        # Check dbNSFP pathogenicity predictions
        pathogenic_indicators = []

        # Check SIFT (low scores = deleterious)
        if "dbNSFP_SIFT_pred" in annotation.dbNSFP_fields:
            sift_pred = annotation.dbNSFP_fields["dbNSFP_SIFT_pred"]
            if "D" in sift_pred:  # Deleterious
                pathogenic_indicators.append(True)

        # Check PolyPhen (high scores = damaging)
        if "dbNSFP_Polyphen2_HDIV_pred" in annotation.dbNSFP_fields:
            poly_pred = annotation.dbNSFP_fields["dbNSFP_Polyphen2_HDIV_pred"]
            if "D" in poly_pred:  # Damaging
                pathogenic_indicators.append(True)

        # Check CADD score (high scores = pathogenic)
        if "dbNSFP_CADD_phred" in annotation.dbNSFP_fields:
            try:
                cadd_score = float(annotation.dbNSFP_fields["dbNSFP_CADD_phred"])
                if cadd_score > 20:  # CADD > 20 suggests pathogenic
                    pathogenic_indicators.append(True)
            except Exception:
                pass

        # Moderate impact with multiple pathogenic indicators
        if annotation.impact == "MODERATE" and len(pathogenic_indicators) >= 2:
            return True

        return False

    def generate_genotype(self, gene: str, annotation: AnnotationTemplate, is_case: bool) -> str:
        """Generate realistic genotype based on gene, annotation, and sample type."""
        genotype = "0/0"  # Default: no variant

        # Determine if this annotation suggests pathogenic variant
        is_pathogenic = self.is_pathogenic_annotation(annotation)

        # Determine variant presence probability
        should_have_variant = False

        if gene in self.DISEASE_GENES:
            if is_case:
                prob = (
                    self.PROBABILITIES["case_disease_pathogenic"]
                    if is_pathogenic
                    else self.PROBABILITIES["case_disease_benign"]
                )
            else:
                prob = (
                    self.PROBABILITIES["control_disease_pathogenic"]
                    if is_pathogenic
                    else self.PROBABILITIES["control_disease_benign"]
                )
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
        pathogenicity = "pathogenic" if is_pathogenic else "benign"
        self.variant_stats[f"{gene}_{sample_type}_{pathogenicity}"][genotype] += 1

        # Return formatted genotype with realistic quality scores
        if genotype == "0/0":
            return f"{genotype}:55:0,55:99"
        elif genotype == "0/1":
            return f"{genotype}:65:28,32:95"
        else:  # 1/1
            return f"{genotype}:45:0,45:99"

    def create_enhanced_vcf_file(self, output_path: Path) -> Dict[str, int]:
        """Create VCF with realistic annotations from sampled data."""
        logger.info(f"Creating enhanced VCF with real annotations: {output_path}")

        stats = {
            "total_variants": 0,
            "disease_variants": 0,
            "control_variants": 0,
            "high_impact_variants": 0,
            "moderate_impact_variants": 0,
        }

        # Generate variants with realistic annotations
        variant_data = self.generate_realistic_variants_with_annotations()

        with open(output_path, "w") as f:
            # Write comprehensive VCF header
            f.write("##fileformat=VCFv4.2\n")
            f.write("##reference=hg19\n")
            f.write(f"##fileDate={datetime.now().strftime('%Y%m%d')}\n")
            f.write("##source=VariantCentrifuge_EnhancedTestDataGenerator_v2.0\n")

            # Add contig lines
            chromosomes = set()
            for gene, (chrom, start, end, strand) in self.GENE_COORDINATES.items():
                chromosomes.add(chrom)

            for chrom in sorted(chromosomes, key=lambda x: int(x) if x.isdigit() else ord(x[0])):
                f.write(f"##contig=<ID={chrom}>\n")

            # INFO field definitions for rich annotations
            info_fields = [
                '##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations: Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID | Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length | AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO">',
                '##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">',
                '##INFO=<ID=SNP,Number=0,Type=Flag,Description="Variant is a SNP">',
                '##INFO=<ID=VARTYPE,Number=1,Type=String,Description="Variant type">',
            ]

            # Dynamically generate annotation field definitions from templates
            annotation_fields = []
            all_field_names = set()

            # Collect all field names from templates
            for templates in self.annotation_sampler.templates.values():
                for template in templates:
                    for key in template.dbNSFP_fields.keys():
                        # Skip problematic keys
                        if "|" in key or "&" in key or len(key) > 50:
                            continue
                        clean_key = key.replace("|", "_").replace("&", "_").replace("=", "_")
                        all_field_names.add(clean_key)
                    for key in template.other_fields.keys():
                        # Skip problematic keys
                        if "|" in key or "&" in key or len(key) > 50:
                            continue
                        clean_key = key.replace("|", "_").replace("&", "_").replace("=", "_")
                        all_field_names.add(clean_key)

            # Add commonly expected fields
            expected_fields = {
                "splice_dbscSNV_rf_score": "Float",
                "splice_dbscSNV_ada_score": "Float",
                "splice_spidex_dpsi_zscore": "Float",
                "dbNSFP_gnomAD_exomes_AC": "Integer",
                "dbNSFP_gnomAD_genomes_AC": "Integer",
                "dbNSFP_ALFA_Total_AC": "Integer",
                "dbNSFP_clinvar_clnsig": "String",
                "ClinVar_CLNSIG": "String",
                "hgmd_CLASS": "String",
            }

            all_field_names.update(expected_fields.keys())

            # Generate field definitions
            for field in sorted(all_field_names):
                if field.startswith("dbNSFP_") and ("AC" in field or "AN" in field):
                    field_type = "Integer"
                elif field.startswith("dbNSFP_") and (
                    "AF" in field or "score" in field or "rankscore" in field
                ):
                    field_type = "Float"
                elif field.startswith("splice_") and "score" in field:
                    field_type = "Float"
                elif field in expected_fields:
                    field_type = expected_fields[field]
                else:
                    field_type = "String"

                annotation_fields.append(
                    f'##INFO=<ID={field},Number=.,Type={field_type},Description="{field} annotation">'
                )

            # Add core dbNSFP fields with proper descriptions
            core_fields = [
                '##INFO=<ID=dbNSFP_CADD_phred,Number=.,Type=Float,Description="CADD phred-like score">',
                '##INFO=<ID=dbNSFP_SIFT_pred,Number=.,Type=String,Description="SIFT prediction">',
                '##INFO=<ID=dbNSFP_Polyphen2_HDIV_pred,Number=.,Type=String,Description="PolyPhen2 HDIV prediction">',
                '##INFO=<ID=dbNSFP_REVEL_score,Number=.,Type=Float,Description="REVEL score">',
                '##INFO=<ID=dbNSFP_GERP___RS,Number=.,Type=Float,Description="GERP++ rejection score">',
            ]

            # Replace any duplicate entries from core_fields
            final_annotation_fields = []
            core_field_ids = set()
            for field in core_fields:
                field_id = field.split("ID=")[1].split(",")[0]
                core_field_ids.add(field_id)
                final_annotation_fields.append(field)

            # Add remaining annotation fields that aren't in core_fields
            for field in annotation_fields:
                field_id = field.split("ID=")[1].split(",")[0]
                if field_id not in core_field_ids:
                    final_annotation_fields.append(field)

            annotation_fields = final_annotation_fields

            for field in info_fields + annotation_fields:
                f.write(field + "\n")

            # FORMAT field definitions
            format_fields = [
                '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">',
                '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read Depth">',
                '##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">',
                '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype Quality">',
            ]

            for field in format_fields:
                f.write(field + "\n")

            # Write sample header
            header_fields = [
                "#CHROM",
                "POS",
                "ID",
                "REF",
                "ALT",
                "QUAL",
                "FILTER",
                "INFO",
                "FORMAT",
            ]
            header = "\t".join(header_fields + self.all_samples) + "\n"
            f.write(header)

            # Write variants with realistic annotations
            for chrom, pos, ref, alt, gene, annotation in variant_data:
                stats["total_variants"] += 1

                if gene in self.DISEASE_GENES:
                    stats["disease_variants"] += 1
                else:
                    stats["control_variants"] += 1

                if annotation.impact == "HIGH":
                    stats["high_impact_variants"] += 1
                elif annotation.impact == "MODERATE":
                    stats["moderate_impact_variants"] += 1

                # Generate genotypes for all samples using realistic annotation-based logic
                genotypes = []
                for sample in self.all_samples:
                    is_case = sample in self.case_samples
                    genotype = self.generate_genotype(gene, annotation, is_case)
                    genotypes.append(genotype)

                # Create realistic INFO field from annotation (including AC from genotypes)
                info = self.create_info_field(annotation, gene, genotypes)

                # Write variant line
                fields = [chrom, pos, ".", ref, alt, "100", "PASS", info, "GT:DP:AD:GQ"]
                line = "\t".join(fields + genotypes) + "\n"
                f.write(line)

        logger.info(f"Created enhanced VCF with {stats['total_variants']} variants")
        logger.info(f"  - Disease variants: {stats['disease_variants']}")
        logger.info(f"  - Control variants: {stats['control_variants']}")
        logger.info(f"  - High impact: {stats['high_impact_variants']}")
        logger.info(f"  - Moderate impact: {stats['moderate_impact_variants']}")

        return stats

    # Include all the helper methods from the original class
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

    def create_phenotype_files(self, output_dir: Path) -> None:
        """Create comprehensive phenotype files for all test scenarios."""
        logger.info("Creating phenotype files")

        # Create all phenotype file variants
        self._create_basic_phenotype_file(output_dir / "phenotypes_basic.csv")
        self._create_extended_phenotype_file(output_dir / "phenotypes_extended.csv")
        self._create_tsv_phenotype_file(output_dir / "phenotypes.tsv")
        self._create_alternative_phenotype_file(output_dir / "phenotypes_alt_columns.csv")
        self._create_incomplete_phenotype_file(output_dir / "phenotypes_incomplete.csv")

    def _create_basic_phenotype_file(self, output_path: Path) -> None:
        """Create basic phenotype file similar to GCKD format."""
        with open(output_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(
                [
                    "MNPPSD",
                    "CGRSequenceID",
                    "SampleID",
                    "Sex",
                    "category_term",
                    "identifier",
                    "name",
                    "link",
                ]
            )

            for i, sample in enumerate(self.all_samples):
                is_case = sample in self.case_samples
                sex = random.choice(["male", "female"])

                hpo_term = (
                    random.choice(self.HPO_TERMS["case_terms"])
                    if is_case
                    else random.choice(self.HPO_TERMS["control_terms"])
                )
                name = self._get_hpo_name(hpo_term)

                writer.writerow(
                    [
                        f"mnp{i:03d}",
                        f"{300000 + i}",
                        sample,
                        sex,
                        "hpo_identifier",
                        hpo_term,
                        name,
                        f"https://hpo.jax.org/app/browse/term/{hpo_term}",
                    ]
                )

    def _create_extended_phenotype_file(self, output_path: Path) -> None:
        """Create phenotype file with multiple HPO terms per sample."""
        with open(output_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(
                [
                    "MNPPSD",
                    "CGRSequenceID",
                    "SampleID",
                    "Sex",
                    "category_term",
                    "identifier",
                    "name",
                    "link",
                ]
            )

            for i, sample in enumerate(self.all_samples):
                is_case = sample in self.case_samples
                sex = random.choice(["male", "female"])
                num_terms = random.randint(1, 3)

                if is_case:
                    available_terms = self.HPO_TERMS["case_terms"] + self.HPO_TERMS["mixed_terms"]
                else:
                    available_terms = (
                        self.HPO_TERMS["control_terms"] + self.HPO_TERMS["mixed_terms"]
                    )

                selected_terms = random.sample(
                    available_terms, min(num_terms, len(available_terms))
                )

                for hpo_term in selected_terms:
                    writer.writerow(
                        [
                            f"mnp{i:03d}",
                            f"{300000 + i}",
                            sample,
                            sex,
                            "hpo_identifier",
                            hpo_term,
                            self._get_hpo_name(hpo_term),
                            f"https://hpo.jax.org/app/browse/term/{hpo_term}",
                        ]
                    )

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
            writer.writerow(
                [
                    "patient_id",
                    "sequencing_id",
                    "sample_name",
                    "gender",
                    "term_type",
                    "hpo_id",
                    "description",
                    "url",
                ]
            )

            for i, sample in enumerate(self.all_samples):
                is_case = sample in self.case_samples
                sex = random.choice(["male", "female"])

                hpo_term = (
                    random.choice(self.HPO_TERMS["case_terms"])
                    if is_case
                    else random.choice(self.HPO_TERMS["control_terms"])
                )

                writer.writerow(
                    [
                        f"patient_{i:03d}",
                        f"{400000 + i}",
                        sample,
                        sex,
                        "hpo_identifier",
                        hpo_term,
                        self._get_hpo_name(hpo_term),
                        f"https://hpo.jax.org/app/browse/term/{hpo_term}",
                    ]
                )

    def _create_incomplete_phenotype_file(self, output_path: Path) -> None:
        """Create phenotype file missing some samples (for error testing)."""
        with open(output_path, "w", newline="") as f:
            writer = csv.writer(f)
            writer.writerow(
                [
                    "MNPPSD",
                    "CGRSequenceID",
                    "SampleID",
                    "Sex",
                    "category_term",
                    "identifier",
                    "name",
                    "link",
                ]
            )

            included_samples = self.all_samples[: int(0.8 * len(self.all_samples))]

            for i, sample in enumerate(included_samples):
                is_case = sample in self.case_samples
                sex = random.choice(["male", "female"])

                hpo_term = (
                    random.choice(self.HPO_TERMS["case_terms"])
                    if is_case
                    else random.choice(self.HPO_TERMS["control_terms"])
                )

                writer.writerow(
                    [
                        f"mnp{i:03d}",
                        f"{300000 + i}",
                        sample,
                        sex,
                        "hpo_identifier",
                        hpo_term,
                        self._get_hpo_name(hpo_term),
                        f"https://hpo.jax.org/app/browse/term/{hpo_term}",
                    ]
                )

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

        (output_dir / "case_hpo_terms.txt").write_text(
            "\n".join(self.HPO_TERMS["case_terms"]) + "\n"
        )
        (output_dir / "control_hpo_terms.txt").write_text(
            "\n".join(self.HPO_TERMS["control_terms"]) + "\n"
        )
        (output_dir / "mixed_hpo_terms.txt").write_text(
            "\n".join(self.HPO_TERMS["mixed_terms"]) + "\n"
        )

    def create_gene_files(self, output_dir: Path) -> None:
        """Create gene list files."""
        logger.info("Creating gene list files")

        (output_dir / "test_genes.txt").write_text("\n".join(sorted(self.ALL_GENES)) + "\n")
        (output_dir / "disease_genes.txt").write_text("\n".join(sorted(self.DISEASE_GENES)) + "\n")
        (output_dir / "control_genes.txt").write_text("\n".join(sorted(self.CONTROL_GENES)) + "\n")

    def write_enhanced_statistics(self, output_dir: Path, vcf_stats: Dict[str, int]) -> None:
        """Write comprehensive statistics about the enhanced dataset."""
        logger.info("Writing enhanced dataset statistics")

        # Calculate detailed genotype distribution
        genotype_summary = {}
        for key, genotype_counts in self.variant_stats.items():
            total = sum(genotype_counts.values())
            if total > 0:
                genotype_summary[key] = {
                    "total_samples": total,
                    "ref_only (0/0)": genotype_counts.get("0/0", 0),
                    "heterozygous (0/1)": genotype_counts.get("0/1", 0),
                    "homozygous (1/1)": genotype_counts.get("1/1", 0),
                    "variant_rate": (genotype_counts.get("0/1", 0) + genotype_counts.get("1/1", 0))
                    / total,
                }

        stats = {
            "generation_info": {
                "date": datetime.now().isoformat(),
                "generator_version": "enhanced_v2.0",
                "description": "Enhanced gene burden test dataset with real annotation sampling",
                "annotation_source": self.annotation_sampler.source_file,
            },
            "sample_composition": {
                "total_samples": len(self.all_samples),
                "case_samples": len(self.case_samples),
                "control_samples": len(self.control_samples),
                "case_sample_ids": self.case_samples,
                "control_sample_ids": self.control_samples,
            },
            "variant_statistics": vcf_stats,
            "annotation_statistics": {
                "total_effect_types": len(self.annotation_sampler.templates),
                "available_effects": list(self.annotation_sampler.templates.keys()),
                "high_impact_effects": self.annotation_sampler.get_high_impact_effects(),
                "moderate_impact_effects": self.annotation_sampler.get_moderate_impact_effects(),
            },
            "genotype_distribution": genotype_summary,
            "test_genes": {
                "disease_genes": list(self.DISEASE_GENES),
                "control_genes": list(self.CONTROL_GENES),
                "total_genes": len(self.ALL_GENES),
            },
            "hpo_terms": self.HPO_TERMS,
            "probabilities": self.PROBABILITIES,
        }

        with open(output_dir / "enhanced_dataset_statistics.json", "w") as f:
            json.dump(stats, f, indent=2)

    def write_enhanced_readme(self, output_dir: Path) -> None:
        """Write comprehensive README for the enhanced dataset."""
        readme_content = """# Enhanced Gene Burden Test Dataset with Real Annotations

Generated on: {datetime.now().isoformat()}

## Overview

This enhanced test dataset provides comprehensive coverage for testing gene burden analysis
functionality in VariantCentrifuge, featuring **realistic annotations sampled from real genomic data**.

### Key Enhancements

- **Real Annotation Sampling**: Authentic SnpEff annotations extracted from {self.annotation_sampler.source_file}
- **Realistic Pathogenicity Prediction**: Uses actual dbNSFP scores (CADD, SIFT, PolyPhen2, REVEL)
- **Diverse Variant Effects**: {len(self.annotation_sampler.templates)} unique effect types from real data
- **Controlled Test Outcomes**: Maintains probabilistic genotype distributions for reliable testing

## Dataset Composition

- **Total Samples**: {len(self.all_samples)} ({len(self.case_samples)} cases, {len(self.control_samples)} controls)
- **Disease Genes**: {', '.join(sorted(self.DISEASE_GENES))}
- **Control Genes**: {', '.join(sorted(self.CONTROL_GENES))}
- **Annotation Source**: Real variants from {self.annotation_sampler.source_file}

## Enhanced Features

### Realistic Annotations
- Authentic ANN fields from SnpEff annotations
- Real dbNSFP pathogenicity predictions
- Actual genomic conservation scores (GERP, PhastCons)
- Population frequency data (gnomAD)

### Effect Type Diversity
Available effect types from real data:
{chr(10).join(f"- {effect}" for effect in sorted(self.annotation_sampler.get_effect_types()))}

### Pathogenicity-Aware Genotype Generation
Genotype probabilities are adjusted based on:
- Variant impact level (HIGH, MODERATE, LOW, MODIFIER)
- CADD phred scores (>20 = likely pathogenic)
- SIFT predictions (D = deleterious)
- PolyPhen2 predictions (D = damaging)

## Files Generated

### Core Data Files
- `enhanced_test_data.vcf.gz` - VCF with realistic annotations
- `enhanced_test_data.vcf.gz.tbi` - Tabix index

### Test Infrastructure
- `run_enhanced_tests.sh` - Test runner with real annotation validation
- `enhanced_dataset_statistics.json` - Detailed statistics including annotation metrics

[Rest of files same as comprehensive dataset...]

## Expected Results

When running gene burden analysis on this enhanced dataset:

- **Disease genes** should show **significant enrichment** in cases
- **High-impact variants** in disease genes should have stronger association
- **Annotation-based filtering** should work with realistic dbNSFP scores
- **Effect prioritization** should properly rank variant consequences

## Validation

The enhanced dataset provides:
- Realistic annotation diversity for testing annotation-based filters
- Authentic pathogenicity scores for testing scoring algorithms
- Controlled genotype distributions ensuring reliable statistical results
- Comprehensive test coverage across all gene burden analysis methods

Check `enhanced_dataset_statistics.json` for detailed metrics on annotation sampling and effect distributions.
"""

        (output_dir / "README.md").write_text(readme_content)

    def generate_enhanced_dataset(self, output_dir: Path) -> None:
        """Generate the complete enhanced test dataset."""
        logger.info(f"Generating enhanced gene burden test dataset in {output_dir}")

        # Create output directory
        output_dir.mkdir(parents=True, exist_ok=True)

        # Generate enhanced components
        vcf_stats = self.create_enhanced_vcf_file(output_dir / "enhanced_test_data.vc")
        self.create_sample_files(output_dir)
        self.create_phenotype_files(output_dir)
        self.create_hpo_term_files(output_dir)
        self.create_gene_files(output_dir)
        self.write_enhanced_statistics(output_dir, vcf_stats)
        self.write_enhanced_readme(output_dir)

        # Compress VCF
        logger.info("Compressing enhanced VCF file...")
        os.system(
            f"cd {output_dir} && bgzip -f enhanced_test_data.vcf && tabix -p vcf enhanced_test_data.vcf.gz"
        )

        logger.info(f"✓ Enhanced test dataset generated in {output_dir}")
        logger.info("✓ Features realistic annotations sampled from real genomic data")
        logger.info("✓ Provides controlled test outcomes with authentic variant effects")


def main():
    """Generate enhanced test dataset with real annotation sampling."""
    parser = argparse.ArgumentParser(
        description="Generate enhanced gene burden test dataset with real annotations",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )

    parser.add_argument(
        "--output-dir",
        default="output",
        help="Output directory for enhanced test dataset (default: output)",
    )
    parser.add_argument(
        "--annotation-source",
        default=None,
        help="Source file with annotations (default: local anonymized templates)",
    )
    parser.add_argument(
        "--seed", type=int, default=42, help="Random seed for reproducibility (default: 42)"
    )

    args = parser.parse_args()

    # Resolve annotation source path (if provided)
    annotation_source = None
    if args.annotation_source:
        annotation_source = Path(args.annotation_source)
        if not annotation_source.is_absolute():
            # Look in current directory and parent directories
            script_dir = Path(__file__).parent
            potential_paths = [
                Path.cwd() / annotation_source,
                script_dir / annotation_source,
                script_dir.parent.parent.parent / annotation_source,
            ]

            for path in potential_paths:
                if path.exists():
                    annotation_source = path
                    break
            else:
                print(f"Error: Could not find annotation source file: {args.annotation_source}")
                print("Searched in:")
                for path in potential_paths:
                    print(f"  - {path}")
                sys.exit(1)

        annotation_source = str(annotation_source)

    # Initialize enhanced generator
    try:
        generator = EnhancedTestDataGenerator(annotation_source, seed=args.seed)
    except Exception as e:
        print(f"Error initializing generator: {e}")
        sys.exit(1)

    # Generate enhanced dataset
    output_dir = Path(args.output_dir)
    generator.generate_enhanced_dataset(output_dir)

    print("\n🎉 Enhanced gene burden test dataset generated!")
    print(f"📁 Location: {output_dir.absolute()}")
    print(f"🧬 Annotations sampled from: {annotation_source or 'local anonymized templates'}")
    print(f"📊 Effect types available: {len(generator.annotation_sampler.templates)}")
    print("📖 See README.md for detailed usage instructions")


if __name__ == "__main__":
    main()
