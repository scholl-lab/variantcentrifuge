#!/usr/bin/env python3
"""
Call variants using GATK HaplotypeCaller on GIAB BAM files.

This script uses GATK to call variants in the gene regions from the extracted
GIAB HG19 exome BAM files using the human_g1k_v37_decoy reference genome.
"""

import argparse
import logging
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

# Configure logging
log_dir = Path("logs")
log_dir.mkdir(exist_ok=True)
log_file = log_dir / f"variant_calling_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.FileHandler(log_file), logging.StreamHandler()],
)
logger = logging.getLogger(__name__)

# Default paths
DEFAULT_DATA_DIR = Path("data")
DEFAULT_REFERENCE_DIR = Path("reference")
DEFAULT_OUTPUT_DIR = Path("variants")
DEFAULT_BED_FILE = Path(__file__).parent / "gene_regions_GRCh37.bed"


def check_dependencies():
    """Check if required tools are available."""
    required_tools = ["gatk", "samtools"]

    for tool in required_tools:
        try:
            result = subprocess.run([tool, "--version"], capture_output=True, check=True)
            if tool == "gatk":
                # GATK version output goes to stderr
                version_output = result.stderr.decode() if result.stderr else result.stdout.decode()
                version_text = version_output.split()[0] if version_output else "Unknown version"
                logger.info(f"Found GATK: {version_text}")
            else:
                version_output = result.stdout.decode().split("\n")[0]
                logger.info(f"Found {tool}: {version_output}")
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.error(f"{tool} not found. Please install {tool} and ensure it's in PATH.")
            return False

    return True


def find_bam_files(data_dir):
    """Find all GIAB BAM files in the data directory."""
    bam_files = []
    for bam_path in data_dir.glob("*.HG19.exome.subset.bam"):
        bai_path = bam_path.with_suffix(".bam.bai")
        if bai_path.exists():
            bam_files.append(bam_path)
        else:
            logger.warning(f"Missing index for {bam_path.name}, skipping")

    if not bam_files:
        logger.error(f"No indexed BAM files found in {data_dir}")
        return None

    logger.info(f"Found {len(bam_files)} indexed BAM files:")
    for bam in sorted(bam_files):
        size_mb = bam.stat().st_size / (1024 * 1024)
        logger.info(f"  - {bam.name} ({size_mb:.1f} MB)")

    return sorted(bam_files)


def find_reference_files(reference_dir):
    """Find reference genome files."""
    fasta_file = reference_dir / "human_g1k_v37_decoy.fasta"

    if not fasta_file.exists():
        logger.error(f"Reference FASTA not found: {fasta_file}")
        logger.info("Run: python download_reference_genome.py")
        return None, None

    # Check for FASTA index
    fai_file = fasta_file.with_suffix(".fasta.fai")
    if not fai_file.exists():
        logger.info(f"FASTA index not found, creating: {fai_file}")
        if not create_fasta_index(fasta_file):
            return None, None

    # Check for sequence dictionary
    dict_file = fasta_file.with_suffix(".dict")
    if not dict_file.exists():
        logger.info(f"Sequence dictionary not found, creating: {dict_file}")
        logger.info(
            "Tip: Use 'python download_reference_genome.py' to download the pre-built dictionary"
        )
        if not create_sequence_dictionary(fasta_file):
            return None, None

    logger.info("Reference files ready:")
    logger.info(f"  - FASTA: {fasta_file}")
    logger.info(f"  - Index: {fai_file}")
    logger.info(f"  - Dictionary: {dict_file}")

    return fasta_file, fai_file


def create_fasta_index(fasta_file):
    """Create FASTA index using samtools."""
    logger.info(f"Creating FASTA index: {fasta_file.name}")

    cmd = ["samtools", "faidx", str(fasta_file)]
    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        logger.error(f"Failed to create FASTA index: {result.stderr}")
        return False

    logger.info("FASTA index created successfully")
    return True


def create_sequence_dictionary(fasta_file):
    """Create sequence dictionary using GATK."""
    dict_file = fasta_file.with_suffix(".dict")
    logger.info(f"Creating sequence dictionary: {dict_file.name}")

    cmd = ["gatk", "CreateSequenceDictionary", "-R", str(fasta_file), "-O", str(dict_file)]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        logger.error(f"Failed to create sequence dictionary: {result.stderr}")
        return False

    logger.info("Sequence dictionary created successfully")
    return True


def sort_bed_file(bed_file, output_dir):
    """Sort BED file by chromosome and position."""
    sorted_bed = output_dir / "gene_regions_sorted.bed"

    if sorted_bed.exists():
        logger.info(f"Sorted BED file already exists: {sorted_bed}")
        return sorted_bed

    logger.info(f"Sorting BED file: {bed_file} -> {sorted_bed}")

    # Read BED file and sort entries
    bed_entries = []
    with open(bed_file) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                parts = line.split("\t")
                if len(parts) >= 3:
                    chrom = parts[0]
                    start = int(parts[1])
                    end = int(parts[2])
                    bed_entries.append((chrom, start, end, line))

    # Sort by chromosome (convert to int if numeric), then by start position
    def sort_key(entry):
        chrom = entry[0]
        try:
            # Handle numeric chromosomes
            if chrom.isdigit():
                return (0, int(chrom), entry[1])
            elif chrom == "X":
                return (1, 23, entry[1])
            elif chrom == "Y":
                return (1, 24, entry[1])
            elif chrom == "MT" or chrom == "M":
                return (1, 25, entry[1])
            else:
                return (2, chrom, entry[1])
        except Exception:
            return (2, chrom, entry[1])

    bed_entries.sort(key=sort_key)

    # Write sorted BED file
    with open(sorted_bed, "w") as f:
        f.write("# Sorted gene regions for GIAB test data - GRCh37/HG19\n")
        f.write("# Format: chr\tstart\tend\tgene_name\n")
        for _, _, _, line in bed_entries:
            f.write(line + "\n")

    logger.info(f"Sorted BED file created with {len(bed_entries)} regions")
    return sorted_bed


def convert_bed_to_interval_list(bed_file, reference_fasta, output_dir):
    """Convert BED file to GATK interval list format."""
    interval_list = output_dir / "gene_regions.interval_list"

    if interval_list.exists():
        logger.info(f"Interval list already exists: {interval_list}")
        return interval_list

    # Sort BED file first
    sorted_bed = sort_bed_file(bed_file, output_dir)

    logger.info(f"Converting sorted BED to interval list: {sorted_bed} -> {interval_list}")

    cmd = [
        "gatk",
        "BedToIntervalList",
        "-I",
        str(sorted_bed),
        "-O",
        str(interval_list),
        "-SD",
        str(reference_fasta.with_suffix(".dict")),
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        logger.error(f"Failed to convert BED to interval list: {result.stderr}")
        return None

    logger.info("Interval list created successfully")
    return interval_list


def call_variants_multisample(
    bam_files, reference_fasta, interval_list, output_dir, java_memory="4g"
):
    """Call variants for all samples together using GATK HaplotypeCaller."""
    output_vcf = output_dir / "HG002-trio_subset_haplotypecaller.vcf.gz"

    # Check if output already exists
    if output_vcf.exists():
        logger.info(f"Output already exists: {output_vcf}")
        return output_vcf

    logger.info(f"Starting multi-sample variant calling for {len(bam_files)} samples")

    # Build command with multiple -I inputs
    cmd = [
        "gatk",
        "--java-options",
        f"-Xmx{java_memory}",
        "HaplotypeCaller",
        "-R",
        str(reference_fasta),
        "-L",
        str(interval_list),
        "-O",
        str(output_vcf),
    ]

    # Add all BAM files as inputs
    for bam_file in bam_files:
        cmd.extend(["-I", str(bam_file)])

    logger.info(f"Running GATK HaplotypeCaller with {len(bam_files)} samples:")
    for bam_file in bam_files:
        logger.info(f"  - {bam_file.name}")
    logger.info(
        f"Command: gatk --java-options -Xmx{java_memory} HaplotypeCaller "
        f"-R [ref] -L [intervals] -O {output_vcf.name} [multiple -I inputs]"
    )

    start_time = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    end_time = time.time()

    if result.returncode != 0:
        logger.error(f"GATK HaplotypeCaller failed: {result.stderr}")
        return None

    # Check output file exists and is not empty
    if not output_vcf.exists():
        logger.error(f"Output VCF not created: {output_vcf}")
        return None

    file_size_mb = output_vcf.stat().st_size / (1024 * 1024)
    duration_minutes = (end_time - start_time) / 60

    logger.info("Multi-sample variant calling completed:")
    logger.info(f"  - Output: {output_vcf.name} ({file_size_mb:.1f} MB)")
    logger.info(f"  - Duration: {duration_minutes:.1f} minutes")

    return output_vcf


def main():
    """Call variants on GIAB BAM files using GATK HaplotypeCaller."""
    parser = argparse.ArgumentParser(
        description="Call variants on GIAB BAM files using GATK HaplotypeCaller"
    )
    parser.add_argument(
        "--data-dir", type=Path, default=DEFAULT_DATA_DIR, help="Directory containing BAM files"
    )
    parser.add_argument(
        "--reference-dir",
        type=Path,
        default=DEFAULT_REFERENCE_DIR,
        help="Directory containing reference genome",
    )
    parser.add_argument(
        "--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR, help="Output directory for VCF files"
    )
    parser.add_argument(
        "--bed-file", type=Path, default=DEFAULT_BED_FILE, help="BED file with regions to call"
    )
    parser.add_argument("--java-memory", default="4g", help="Java memory allocation for GATK")

    args = parser.parse_args()

    logger.info(f"Starting GATK variant calling - Log file: {log_file}")
    logger.info(
        f"Parameters: data_dir={args.data_dir}, reference_dir={args.reference_dir}, "
        f"output_dir={args.output_dir}"
    )

    # Check dependencies
    if not check_dependencies():
        sys.exit(1)

    # Find input files
    bam_files = find_bam_files(args.data_dir)
    if not bam_files:
        sys.exit(1)

    reference_fasta, reference_fai = find_reference_files(args.reference_dir)
    if not reference_fasta:
        sys.exit(1)

    if not args.bed_file.exists():
        logger.error(f"BED file not found: {args.bed_file}")
        sys.exit(1)

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Convert BED to interval list
    interval_list = convert_bed_to_interval_list(args.bed_file, reference_fasta, args.output_dir)
    if not interval_list:
        sys.exit(1)

    # Call variants for all samples together
    output_vcf = call_variants_multisample(
        bam_files, reference_fasta, interval_list, args.output_dir, args.java_memory
    )

    # Final summary
    logger.info("\n" + "=" * 60)
    logger.info("VARIANT CALLING SUMMARY")
    logger.info("=" * 60)

    if output_vcf:
        size_mb = output_vcf.stat().st_size / (1024 * 1024)
        logger.info(f"✓ Multi-sample VCF created: {output_vcf.name} ({size_mb:.1f} MB)")
        logger.info(f"Samples processed: {len(bam_files)}")
    else:
        logger.error("✗ Variant calling failed")
        sys.exit(1)

    logger.info(f"\nOutput directory: {args.output_dir}")
    logger.info(f"Log saved to: {log_file}")
    logger.info("Done!")


if __name__ == "__main__":
    main()
