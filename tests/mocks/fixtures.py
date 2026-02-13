"""Test fixtures and factory functions."""

import argparse
import tempfile
from pathlib import Path

from variantcentrifuge.pipeline_core import PipelineContext, Workspace


def create_test_context(
    vcf_file: str = "test.vcf",
    output_dir: str = None,
    gene_name: str = "BRCA1",
    config_overrides: dict | None = None,
    **kwargs,
) -> PipelineContext:
    """Create a test PipelineContext with sensible defaults.

    Parameters
    ----------
    vcf_file : str
        Path to VCF file
    output_dir : str, optional
        Output directory (temp dir if not specified)
    gene_name : str
        Gene name for analysis
    config_overrides : dict, optional
        Config values to override
    **kwargs
        Additional arguments for argparse.Namespace

    Returns
    -------
    PipelineContext
        Configured test context
    """
    # Create args with defaults
    default_args = {
        "vcf_file": vcf_file,
        "output_dir": output_dir or tempfile.mkdtemp(),
        "gene_name": gene_name,
        "gene_file": None,
        "output_file": None,
        "config": None,
        "phenotype_file": None,
        "phenotype_sample_column": None,
        "phenotype_value_column": None,
        "ped_file": None,
        "scoring_config_path": None,
        "threads": 1,
    }

    # Update with any kwargs, avoiding duplicates
    default_args.update(kwargs)

    # Create namespace
    args = argparse.Namespace(**default_args)

    # Create workspace
    workspace = Workspace(Path(args.output_dir), "test_run")

    # Create config
    config = {
        "vcf_file": vcf_file,
        "output_dir": args.output_dir,
        "gene_name": gene_name,
        "threads": 1,
        "pipeline_version": "1.0.0-test",
    }

    if config_overrides:
        config.update(config_overrides)

    # Create context
    context = PipelineContext(args=args, config=config, workspace=workspace)

    return context


def create_test_vcf(
    output_path: Path, samples: list[str] = None, variants: list[dict] = None
) -> Path:
    """Create a test VCF file.

    Parameters
    ----------
    output_path : Path
        Where to write the VCF
    samples : list of str, optional
        Sample names (default: Sample1, Sample2, Sample3)
    variants : list of dict, optional
        Variant records

    Returns
    -------
    Path
        Path to created VCF file
    """
    if samples is None:
        samples = ["Sample1", "Sample2", "Sample3"]

    if variants is None:
        variants = [
            {
                "chrom": "chr1",
                "pos": 100,
                "ref": "A",
                "alt": "T",
                "qual": 100,
                "filter": "PASS",
                "info": "AF=0.01",
                "genotypes": ["0/1", "0/0", "0/1"],
            },
            {
                "chrom": "chr1",
                "pos": 200,
                "ref": "G",
                "alt": "C",
                "qual": 100,
                "filter": "PASS",
                "info": "AF=0.02",
                "genotypes": ["0/0", "0/1", "1/1"],
            },
        ]

    # Write VCF
    with open(output_path, "w") as f:
        # Header
        f.write("##fileformat=VCFv4.2\n")
        f.write("##contig=<ID=chr1,length=249250621>\n")
        f.write('##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">\n')
        f.write('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n')

        # Column headers
        cols = ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"]
        cols.extend(samples)
        f.write("\t".join(cols) + "\n")

        # Variants
        for var in variants:
            fields = [
                var["chrom"],
                str(var["pos"]),
                ".",
                var["ref"],
                var["alt"],
                str(var["qual"]),
                var["filter"],
                var["info"],
                "GT",
            ]
            fields.extend(var["genotypes"])
            f.write("\t".join(fields) + "\n")

    return output_path


def create_test_bed(output_path: Path, regions: list[dict] = None) -> Path:
    """Create a test BED file.

    Parameters
    ----------
    output_path : Path
        Where to write the BED
    regions : list of dict, optional
        Region records with chrom, start, end, name

    Returns
    -------
    Path
        Path to created BED file
    """
    if regions is None:
        regions = [
            {"chrom": "chr1", "start": 1000, "end": 2000, "name": "GENE1"},
            {"chrom": "chr1", "start": 5000, "end": 6000, "name": "GENE2"},
            {"chrom": "chr2", "start": 1000, "end": 2000, "name": "GENE3"},
        ]

    with open(output_path, "w") as f:
        for region in regions:
            fields = [
                region["chrom"],
                str(region["start"]),
                str(region["end"]),
                region.get("name", "."),
                region.get("score", "0"),
                region.get("strand", "+"),
            ]
            f.write("\t".join(fields) + "\n")

    return output_path


def create_test_phenotype_file(output_path: Path, phenotypes: dict[str, str] = None) -> Path:
    """Create a test phenotype file.

    Parameters
    ----------
    output_path : Path
        Where to write the file
    phenotypes : dict, optional
        Sample to phenotype mapping

    Returns
    -------
    Path
        Path to created file
    """
    if phenotypes is None:
        phenotypes = {
            "Sample1": "Affected",
            "Sample2": "Unaffected",
            "Sample3": "Affected",
        }

    with open(output_path, "w") as f:
        f.write("Sample\tPhenotype\n")
        for sample, pheno in phenotypes.items():
            f.write(f"{sample}\t{pheno}\n")

    return output_path


def create_test_pedigree_file(output_path: Path, families: list[dict] = None) -> Path:
    """Create a test PED file.

    Parameters
    ----------
    output_path : Path
        Where to write the file
    families : list of dict, optional
        Family information

    Returns
    -------
    Path
        Path to created file
    """
    if families is None:
        families = [
            {
                "family": "FAM001",
                "individual": "Sample1",
                "father": "0",
                "mother": "0",
                "sex": "1",
                "phenotype": "2",
            },
            {
                "family": "FAM001",
                "individual": "Sample2",
                "father": "0",
                "mother": "0",
                "sex": "2",
                "phenotype": "1",
            },
            {
                "family": "FAM001",
                "individual": "Sample3",
                "father": "Sample1",
                "mother": "Sample2",
                "sex": "1",
                "phenotype": "2",
            },
        ]

    with open(output_path, "w") as f:
        for fam in families:
            fields = [
                fam["family"],
                fam["individual"],
                fam["father"],
                fam["mother"],
                fam["sex"],
                fam["phenotype"],
            ]
            f.write("\t".join(fields) + "\n")

    return output_path
