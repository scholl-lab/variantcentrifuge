#!/usr/bin/env python3
"""Download GIAB BAM files and extract specific genomic regions defined in BED files."""

import sys
import subprocess
import argparse
import logging
import hashlib
import json
from pathlib import Path
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed

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

# GIAB BAM URLs and MD5 checksums
GIAB_BAMS = {
    "HG002": {
        "GRCh38": {
            "bam": (
                "ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/"
                "HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/"
                "NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x.bam"
            ),
            "bai": (
                "ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/"
                "HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/"
                "NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.300x.bam.bai"
            ),
            "bam_md5": "a2a85243fba525d485a6ba89d73b46a9",
            "bai_md5": "491d644a55382424c1824e3b664827da",
        },
        "GRCh37": {
            "bam": (
                "ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/"
                "HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/"
                "NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.300x.bam"
            ),
            "bai": (
                "ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/"
                "HG002_NA24385_son/NIST_HiSeq_HG002_Homogeneity-10953946/"
                "NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.300x.bam.bai"
            ),
            "bam_md5": "306da46fc735f3b28ff0d619d6337de2",
            "bai_md5": "98a533624e2070db4a4ad0d6ee4b295e",
        },
    },
    "HG003": {
        "GRCh38": {
            "bam": (
                "ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/"
                "HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/"
                "NHGRI_Illumina300X_AJtrio_novoalign_bams/HG003.GRCh38.300x.bam"
            ),
            "bai": (
                "ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/"
                "HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/"
                "NHGRI_Illumina300X_AJtrio_novoalign_bams/HG003.GRCh38.300x.bam.bai"
            ),
            "bam_md5": "3642fca7bd403073ce7c0b62a5686f66",
            "bai_md5": "e59006c82ad053ca5e40357e26ee2342",
        },
        "GRCh37": {
            "bam": (
                "ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/"
                "HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/"
                "NHGRI_Illumina300X_AJtrio_novoalign_bams/HG003.hs37d5.300x.bam"
            ),
            "bai": (
                "ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/"
                "HG003_NA24149_father/NIST_HiSeq_HG003_Homogeneity-12389378/"
                "NHGRI_Illumina300X_AJtrio_novoalign_bams/HG003.hs37d5.300x.bam.bai"
            ),
            "bam_md5": "478322e9a29c3323f056cdbf8b20ab2f",
            "bai_md5": "4cc888558e946058f4fcaef2086ddd56",
        },
    },
    "HG004": {
        "GRCh38": {
            "bam": (
                "ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/"
                "HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/"
                "NHGRI_Illumina300X_AJtrio_novoalign_bams/HG004.GRCh38.300x.bam"
            ),
            "bai": (
                "ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/"
                "HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/"
                "NHGRI_Illumina300X_AJtrio_novoalign_bams/HG004.GRCh38.300x.bam.bai"
            ),
            "bam_md5": "b0a04174824f9ffa66ac06a16eeb8c73",
            "bai_md5": "0daf310b2d3767f98401b0395583988d",
        },
        "GRCh37": {
            "bam": (
                "ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/"
                "HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/"
                "NHGRI_Illumina300X_AJtrio_novoalign_bams/HG004.hs37d5.300x.bam"
            ),
            "bai": (
                "ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/"
                "HG004_NA24143_mother/NIST_HiSeq_HG004_Homogeneity-14572558/"
                "NHGRI_Illumina300X_AJtrio_novoalign_bams/HG004.hs37d5.300x.bam.bai"
            ),
            "bam_md5": "7b8d1b5a39734633adf0d901e3231cb4",
            "bai_md5": "3e73c04bd7daeebcc8193e3a77d9082a",
        },
    },
}


def calculate_md5(filepath, chunk_size=8192):
    """Calculate MD5 checksum of a file."""
    md5_hash = hashlib.md5()
    with open(filepath, "rb") as f:
        while chunk := f.read(chunk_size):
            md5_hash.update(chunk)
    return md5_hash.hexdigest()


def verify_file_md5(filepath, expected_md5):
    """Verify file MD5 checksum."""
    if not filepath.exists():
        return False

    logger.info(f"Calculating MD5 for {filepath.name}...")
    actual_md5 = calculate_md5(filepath)

    if actual_md5 == expected_md5:
        logger.info(f"MD5 verified for {filepath.name}: {actual_md5}")
        return True
    else:
        logger.error(f"MD5 mismatch for {filepath.name}: expected {expected_md5}, got {actual_md5}")
        return False


def load_download_status():
    """Load download status from JSON file."""
    status_file = Path("download_status.json")
    if status_file.exists():
        with open(status_file, "r") as f:
            return json.load(f)
    return {}


def save_download_status(status):
    """Save download status to JSON file."""
    status_file = Path("download_status.json")
    with open(status_file, "w") as f:
        json.dump(status, f, indent=2)


def download_file(url, output_path, expected_md5=None):
    """Download a file using wget with resume support."""
    logger.info(f"Downloading {output_path.name} from {url}")

    # Check if file already exists and is valid
    if output_path.exists() and expected_md5:
        if verify_file_md5(output_path, expected_md5):
            logger.info(f"{output_path.name} already exists with correct MD5, skipping download")
            return True
        else:
            logger.warning(f"{output_path.name} exists but MD5 mismatch, re-downloading")
            output_path.unlink()

    # Download with wget (resume support)
    cmd = ["wget", "-c", "-O", str(output_path), url]
    result = subprocess.run(cmd, capture_output=True)

    if result.returncode != 0:
        logger.error(f"Failed to download {url}: {result.stderr.decode()}")
        return False

    # Verify MD5 if provided
    if expected_md5:
        if not verify_file_md5(output_path, expected_md5):
            output_path.unlink()  # Remove corrupted file
            return False

    return True


def download_sample_files(sample, assembly, temp_dir):
    """Download BAM and BAI files for a sample."""
    urls = GIAB_BAMS[sample][assembly]
    bam_path = temp_dir / f"{sample}.{assembly}.full.bam"
    bai_path = temp_dir / f"{sample}.{assembly}.full.bam.bai"

    # Load download status
    status = load_download_status()
    sample_key = f"{sample}_{assembly}"

    # Check if already completed
    if sample_key in status and status[sample_key].get("completed", False):
        if bam_path.exists() and bai_path.exists():
            logger.info(f"{sample} {assembly} already downloaded and verified")
            return (sample, bam_path)

    # Download both files
    success = True

    # Download BAI first (smaller file)
    if not download_file(urls["bai"], bai_path, urls.get("bai_md5")):
        success = False

    # Download BAM
    if success and not download_file(urls["bam"], bam_path, urls.get("bam_md5")):
        success = False

    # Update status
    if success:
        status[sample_key] = {
            "completed": True,
            "timestamp": datetime.now().isoformat(),
            "bam_md5": calculate_md5(bam_path),
            "bai_md5": calculate_md5(bai_path),
            "bam_size": bam_path.stat().st_size,
            "bai_size": bai_path.stat().st_size,
        }
        save_download_status(status)
        logger.info(f"Successfully downloaded {sample} {assembly}")

    return (sample, bam_path if success else None)


def extract_regions_from_bam(bam_path, bed_file, output_bam, threads=4):
    """Extract regions from a BAM file using a BED file."""
    # Check if output already exists
    if output_bam.exists() and output_bam.with_suffix(".bam.bai").exists():
        logger.info(f"{output_bam} already exists, checking integrity...")
        # Quick check if BAM is valid
        result = subprocess.run(["samtools", "quickcheck", str(output_bam)], capture_output=True)
        if result.returncode == 0:
            size_mb = output_bam.stat().st_size / (1024 * 1024)
            logger.info(f"{output_bam} is valid ({size_mb:.1f} MB), skipping extraction")
            return True
        else:
            logger.warning(f"{output_bam} is corrupted, re-extracting")
            output_bam.unlink(missing_ok=True)
            output_bam.with_suffix(".bam.bai").unlink(missing_ok=True)

    # Read regions from BED
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

    logger.info(f"Extracting {len(regions)} regions from {bam_path.name}")

    # Create a BED file for samtools view -L
    temp_bed = output_bam.parent / "temp_regions.bed"
    with open(bed_file, "r") as infile, open(temp_bed, "w") as outfile:
        for line in infile:
            if line.strip() and not line.startswith("#"):
                outfile.write(line)

    # Extract regions using BED file
    cmd = [
        "samtools",
        "view",
        "-b",
        "-@",
        str(threads),
        "-L",
        str(temp_bed),
        "-o",
        str(output_bam),
        str(bam_path),
    ]

    result = subprocess.run(cmd, capture_output=True)
    if result.returncode != 0:
        logger.error(f"Failed to extract regions: {result.stderr.decode()}")
        temp_bed.unlink()
        return False

    # Index the output
    logger.info(f"Indexing {output_bam.name}")
    subprocess.run(["samtools", "index", str(output_bam)], check=True)

    # Calculate MD5 of extracted BAM
    output_md5 = calculate_md5(output_bam)
    logger.info(f"MD5 of extracted BAM: {output_md5}")

    # Save extraction info
    extraction_info = {
        "source_bam": str(bam_path),
        "bed_file": str(bed_file),
        "output_bam": str(output_bam),
        "output_md5": output_md5,
        "timestamp": datetime.now().isoformat(),
        "regions": len(regions),
    }

    info_file = output_bam.with_suffix(".extraction_info.json")
    with open(info_file, "w") as f:
        json.dump(extraction_info, f, indent=2)

    # Cleanup
    temp_bed.unlink()

    # Report size
    size_mb = output_bam.stat().st_size / (1024 * 1024)
    logger.info(f"Created {output_bam} ({size_mb:.1f} MB)")

    return True


def main():
    """Download GIAB BAM files and extract gene regions."""
    parser = argparse.ArgumentParser(description="Download GIAB BAM files and extract gene regions")
    parser.add_argument("--output-dir", type=Path, default=Path("data"), help="Output directory")
    parser.add_argument(
        "--samples",
        nargs="+",
        choices=["HG002", "HG003", "HG004", "all"],
        default=["all"],
        help="Samples to download",
    )
    parser.add_argument(
        "--assembly", choices=["GRCh37", "GRCh38"], required=True, help="Genome assembly"
    )
    parser.add_argument("--threads", type=int, default=4, help="Number of threads")
    parser.add_argument(
        "--keep-full", action="store_true", help="Keep full BAM files after extraction"
    )
    parser.add_argument(
        "--skip-md5", action="store_true", help="Skip MD5 verification for faster processing"
    )

    args = parser.parse_args()

    logger.info(f"Starting GIAB download script - Log file: {log_file}")
    logger.info(
        f"Parameters: assembly={args.assembly}, samples={args.samples}, " f"threads={args.threads}"
    )

    # Check dependencies
    for tool in ["samtools", "wget"]:
        try:
            subprocess.run([tool, "--version"], capture_output=True, check=True)
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.error(f"{tool} not found. Please install {tool} first.")
            sys.exit(1)

    # Get BED file
    bed_file = Path(__file__).parent / f"gene_regions_{args.assembly}.bed"
    if not bed_file.exists():
        logger.error(f"BED file not found: {bed_file}")
        sys.exit(1)

    # Get samples
    samples = ["HG002", "HG003", "HG004"] if "all" in args.samples else args.samples

    # Create directories
    args.output_dir.mkdir(parents=True, exist_ok=True)
    temp_dir = args.output_dir / "temp"
    temp_dir.mkdir(exist_ok=True)

    # Download files in parallel
    logger.info(f"Processing {len(samples)} samples...")
    download_tasks = []

    with ThreadPoolExecutor(max_workers=min(len(samples), args.threads)) as executor:
        # Submit download tasks
        for sample in samples:
            output_bam = args.output_dir / f"{sample}.{args.assembly}.bam"

            # Check if final output already exists
            if output_bam.exists() and output_bam.with_suffix(".bam.bai").exists():
                result = subprocess.run(
                    ["samtools", "quickcheck", str(output_bam)], capture_output=True
                )
                if result.returncode == 0:
                    logger.info(f"{output_bam} already exists and is valid, skipping")
                    continue

            future = executor.submit(download_sample_files, sample, args.assembly, temp_dir)
            download_tasks.append(future)

        # Process completed downloads
        for future in as_completed(download_tasks):
            sample, full_bam_path = future.result()
            if full_bam_path:
                logger.info(f"Download complete for {sample}, starting extraction")

                # Extract regions
                output_bam = args.output_dir / f"{sample}.{args.assembly}.bam"
                success = extract_regions_from_bam(
                    full_bam_path, bed_file, output_bam, args.threads
                )

                if success and not args.keep_full:
                    # Remove full BAM files
                    logger.info(f"Removing full BAM files for {sample}")
                    full_bam_path.unlink()
                    full_bam_path.with_suffix(".bam.bai").unlink()
            else:
                logger.error(f"Download failed for {sample}")

    # Cleanup temp directory if empty
    if not args.keep_full and not any(temp_dir.iterdir()):
        temp_dir.rmdir()

    # Final summary
    logger.info("\n" + "=" * 60)
    logger.info("DOWNLOAD SUMMARY")
    logger.info("=" * 60)

    # Report successful extractions
    extracted_bams = list(args.output_dir.glob("*.bam"))
    extracted_bams = [b for b in extracted_bams if "full" not in b.name]

    if extracted_bams:
        logger.info(f"Successfully extracted {len(extracted_bams)} BAM files:")
        total_size = 0
        for bam in sorted(extracted_bams):
            size_mb = bam.stat().st_size / (1024 * 1024)
            total_size += size_mb
            logger.info(f"  - {bam.name}: {size_mb:.1f} MB")
        logger.info(f"Total size: {total_size:.1f} MB")
    else:
        logger.warning("No BAM files were successfully extracted")

    logger.info(f"\nLog saved to: {log_file}")
    logger.info("Done!")


if __name__ == "__main__":
    main()
