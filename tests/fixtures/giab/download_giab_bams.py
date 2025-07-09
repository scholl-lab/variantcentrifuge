#!/usr/bin/env python3
"""Extract specific genomic regions from remote GIAB BAM files using samtools."""

import argparse
import hashlib
import json
import logging
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path

# Configure logging
log_dir = Path("logs")
log_dir.mkdir(exist_ok=True)
log_file = log_dir / f"download_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

# Set up both file and console logging
logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.FileHandler(log_file), logging.StreamHandler()],
)
logger = logging.getLogger(__name__)

# GIAB HG19 Exome BAM URLs
GIAB_EXOME_BAMS = {
    "HG002": {
        "bam": "ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bam",
        "bai": "ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG002-EEogPU_v02-KIT-Av5_AGATGTAC_L008.posiSrt.markDup.bai",
        "bam_md5": "c80f0cab24bfaa504393457b8f7191fa",
        "bai_md5": "d4fea426c3e2e9a71bb92e6526b4df6",
    },
    "HG003": {
        "bam": "ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG003_NA24149_father/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG003-EEogPU_v02-KIT-Av5_TCTTCACA_L008.posiSrt.markDup.bam",
        "bai": "ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG003_NA24149_father/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG003-EEogPU_v02-KIT-Av5_TCTTCACA_L008.posiSrt.markDup.bai",
        "bam_md5": "79ce0edc9363710030e973bd12af5423",
        "bai_md5": "0512b54f2b66f3a04653be9f5975ae7",
    },
    "HG004": {
        "bam": "ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG004_NA24143_mother/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG004-EEogPU_v02-KIT-Av5_CCGAAGTA_L008.posiSrt.markDup.bam",
        "bai": "ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG004_NA24143_mother/OsloUniversityHospital_Exome/151002_7001448_0359_AC7F6GANXX_Sample_HG004-EEogPU_v02-KIT-Av5_CCGAAGTA_L008.posiSrt.markDup.bai",
        "bam_md5": "39e90d925fff5ea7dd6dab255f299581",
        "bai_md5": "8914bfb6fa6bd304192f2c9e13e903f4",
    },
}


def calculate_md5(filepath, chunk_size=8192):
    """Calculate MD5 checksum of a file."""
    md5_hash = hashlib.md5()
    with open(filepath, "rb") as f:
        while chunk := f.read(chunk_size):
            md5_hash.update(chunk)
    return md5_hash.hexdigest()


def load_extraction_status():
    """Load extraction status from JSON file."""
    status_file = Path("extraction_status.json")
    if status_file.exists():
        with open(status_file, "r") as f:
            return json.load(f)
    return {}


def save_extraction_status(status):
    """Save extraction status to JSON file."""
    status_file = Path("extraction_status.json")
    with open(status_file, "w") as f:
        json.dump(status, f, indent=2)


def extract_regions_from_remote_bam(sample, bam_url, bed_file, output_bam, threads=4):
    """Extract regions directly from remote BAM URL using samtools."""
    logger.info(f"Starting extraction: {sample} -> {output_bam.name}")

    # Check if output already exists and is valid
    if output_bam.exists() and output_bam.with_suffix(".bam.bai").exists():
        logger.info(f"{output_bam} already exists, checking integrity...")
        result = subprocess.run(["samtools", "quickcheck", str(output_bam)], capture_output=True)
        if result.returncode == 0:
            size_mb = output_bam.stat().st_size / (1024 * 1024)
            logger.info(f"{output_bam} is valid ({size_mb:.1f} MB), skipping extraction")
            return True
        else:
            logger.warning(f"{output_bam} is corrupted, re-extracting")
            output_bam.unlink(missing_ok=True)
            output_bam.with_suffix(".bam.bai").unlink(missing_ok=True)

    # Count regions for logging
    regions = []
    with open(bed_file, "r") as f:
        for line in f:
            if line.strip() and not line.startswith("#"):
                parts = line.strip().split("\t")
                if len(parts) >= 3:
                    regions.append(f"{parts[0]}:{parts[1]}-{parts[2]}")

    if not regions:
        logger.error(f"No regions found in {bed_file}")
        return False

    logger.info(
        f"Extracting {len(regions)} regions from remote BAM (this streams data, no local download)"
    )

    # Build command with individual regions for remote BAM
    cmd = ["samtools", "view", "-b", "-@", str(threads), "-o", str(output_bam), bam_url]

    # Add each region individually
    for region in regions:
        cmd.append(region)

    logger.info(
        f"Running: samtools view -b -@ {threads} -o {output_bam.name} {bam_url} {' '.join(regions)}"
    )
    logger.info("Streaming extraction from remote BAM - this may take several minutes...")

    start_time = time.time()
    result = subprocess.run(cmd, capture_output=True, text=True)
    end_time = time.time()

    if result.returncode != 0:
        logger.error(f"Failed to extract regions from remote BAM: {result.stderr}")
        return False

    logger.info(f"Remote extraction completed successfully in {end_time - start_time:.1f} seconds")

    # Sort the BAM file (required because multiple regions create unsorted output)
    sorted_bam = output_bam.with_suffix(".sorted.bam")
    logger.info(f"Sorting {output_bam.name}...")
    sort_start = time.time()
    sort_result = subprocess.run(
        ["samtools", "sort", "-@", str(threads), "-o", str(sorted_bam), str(output_bam)],
        capture_output=True,
        text=True,
    )
    sort_end = time.time()

    if sort_result.returncode != 0:
        logger.error(f"Failed to sort BAM: {sort_result.stderr}")
        return False

    # Replace unsorted with sorted
    output_bam.unlink()
    sorted_bam.rename(output_bam)
    logger.info(f"Sorting completed successfully in {sort_end - sort_start:.1f} seconds")

    # Index the sorted output
    logger.info(f"Indexing {output_bam.name}")
    index_result = subprocess.run(
        ["samtools", "index", str(output_bam)], capture_output=True, text=True
    )
    if index_result.returncode != 0:
        logger.error(f"Failed to index BAM: {index_result.stderr}")
        return False
    logger.info("Indexing completed successfully")

    # Calculate MD5 of extracted BAM
    logger.info(f"Calculating MD5 for {output_bam.name}...")
    output_md5 = calculate_md5(output_bam)
    logger.info(f"MD5 of extracted BAM: {output_md5}")

    # Save extraction info
    extraction_info = {
        "sample": sample,
        "assembly": "HG19",
        "source_url": bam_url,
        "bed_file": str(bed_file),
        "output_bam": str(output_bam),
        "output_md5": output_md5,
        "timestamp": datetime.now().isoformat(),
        "regions": len(regions),
        "extraction_time_seconds": end_time - start_time,
        "sort_time_seconds": sort_end - sort_start,
    }

    info_file = output_bam.with_suffix(".extraction_info.json")
    with open(info_file, "w") as f:
        json.dump(extraction_info, f, indent=2)

    # Report size
    size_mb = output_bam.stat().st_size / (1024 * 1024)
    logger.info(f"Created {output_bam} ({size_mb:.1f} MB)")

    # Clean up any temporary files that might have been created during processing
    temp_files = [
        Path.cwd() / f"{sample}.bai",
        Path.cwd() / f"{sample}.HG19.exome.bai",
        Path.cwd() / "temp_regions.bed",
    ]

    for temp_file in temp_files:
        if temp_file.exists():
            logger.info(f"Cleaning up temporary file: {temp_file}")
            temp_file.unlink()

    return True


def main():
    """Extract gene regions from remote GIAB HG19 exome BAM files using samtools."""
    parser = argparse.ArgumentParser(
        description="Extract gene regions from remote GIAB HG19 exome BAM files"
    )
    parser.add_argument("--output-dir", type=Path, default=Path("data"), help="Output directory")
    parser.add_argument(
        "--samples",
        nargs="+",
        choices=["HG002", "HG003", "HG004", "all"],
        default=["all"],
        help="Samples to extract",
    )
    parser.add_argument("--threads", type=int, default=4, help="Number of threads for samtools")

    args = parser.parse_args()

    logger.info(f"Starting GIAB HG19 exome extraction script - Log file: {log_file}")
    logger.info(f"Parameters: samples={args.samples}, threads={args.threads}")

    # Check dependencies - only need samtools now
    try:
        subprocess.run(["samtools", "--version"], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        logger.error("samtools not found. Please install samtools first.")
        sys.exit(1)

    # Get BED file - using HG19/GRCh37 coordinates
    bed_file = Path(__file__).parent / "gene_regions_GRCh37.bed"
    if not bed_file.exists():
        logger.error(f"BED file not found: {bed_file}")
        sys.exit(1)

    # Get samples
    samples = ["HG002", "HG003", "HG004"] if "all" in args.samples else args.samples

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Extract regions from remote BAMs in parallel
    logger.info(f"Extracting regions from {len(samples)} remote HG19 exome BAM files...")
    extraction_tasks = []

    with ThreadPoolExecutor(max_workers=min(len(samples), args.threads)) as executor:
        # Submit extraction tasks
        for sample in samples:
            output_bam = args.output_dir / f"{sample}.HG19.exome.subset.bam"

            # Check if final output already exists and is valid
            if output_bam.exists() and output_bam.with_suffix(".bam.bai").exists():
                result = subprocess.run(
                    ["samtools", "quickcheck", str(output_bam)], capture_output=True
                )
                if result.returncode == 0:
                    logger.info(f"{output_bam} already exists and is valid, skipping")
                    continue

            # Get remote BAM URL
            bam_url = GIAB_EXOME_BAMS[sample]["bam"]

            # Submit extraction task
            future = executor.submit(
                extract_regions_from_remote_bam, sample, bam_url, bed_file, output_bam, args.threads
            )
            extraction_tasks.append((future, sample))

        # Process completed extractions
        for future, sample in extraction_tasks:
            success = future.result()
            if success:
                logger.info(f"Extraction complete for {sample}")
            else:
                logger.error(f"Extraction failed for {sample}")

    # Final summary
    logger.info("\n" + "=" * 60)
    logger.info("EXTRACTION SUMMARY")
    logger.info("=" * 60)

    # Report successful extractions
    extracted_bams = list(args.output_dir.glob("*.bam"))

    if extracted_bams:
        logger.info(f"Successfully extracted {len(extracted_bams)} BAM files:")
        total_size = 0
        for bam in sorted(extracted_bams):
            size_mb = bam.stat().st_size / (1024 * 1024)
            total_size += size_mb
            logger.info(f"  - {bam.name}: {size_mb:.1f} MB")
        logger.info(f"Total size: {total_size:.1f} MB")

        # Show total extraction info files
        info_files = list(args.output_dir.glob("*.extraction_info.json"))
        if info_files:
            logger.info(
                f"Extraction metadata saved in {len(info_files)} .extraction_info.json files"
            )
    else:
        logger.warning("No BAM files were successfully extracted")

    # Final cleanup of any stray temporary files in working directory
    logger.info("Performing final cleanup...")
    cleanup_patterns = ["*.bai", "temp_*.bed", "*temp*.bam"]
    cleanup_count = 0

    for pattern in cleanup_patterns:
        for temp_file in Path.cwd().glob(pattern):
            if temp_file.is_file() and temp_file.name not in [f.name for f in extracted_bams]:
                logger.info(f"Removing temporary file: {temp_file}")
                temp_file.unlink()
                cleanup_count += 1

    if cleanup_count > 0:
        logger.info(f"Cleaned up {cleanup_count} temporary files")
    else:
        logger.info("No temporary files to clean up")

    logger.info(f"\nLog saved to: {log_file}")
    logger.info("Done!")


if __name__ == "__main__":
    main()
