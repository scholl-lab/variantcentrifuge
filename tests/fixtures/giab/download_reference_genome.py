#!/usr/bin/env python3
"""
Download human_g1k_v37_decoy reference genome from Broad Institute FTP server.

This script downloads the standard reference genome used in variant calling pipelines:
- human_g1k_v37_decoy.fasta.gz (reference sequence)
- human_g1k_v37_decoy.fasta.fai.gz (FASTA index)

The reference includes decoy sequences to improve alignment accuracy.
"""

import argparse
import gzip
import hashlib
import logging
import shutil
import subprocess
import sys
from datetime import datetime
from pathlib import Path

# Configure logging
log_dir = Path("logs")
log_dir.mkdir(exist_ok=True)
log_file = log_dir / f"reference_download_{datetime.now().strftime('%Y%m%d_%H%M%S')}.log"

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
    handlers=[logging.FileHandler(log_file), logging.StreamHandler()],
)
logger = logging.getLogger(__name__)

# Broad Institute FTP server details
FTP_SERVER = "ftp.broadinstitute.org"
FTP_USERNAME = "gsapubftp-anonymous"
FTP_PASSWORD = ""
FTP_BASE_PATH = "/bundle/b37"

# Files to download
REFERENCE_FILES = {
    "fasta": {
        "remote": "human_g1k_v37_decoy.fasta.gz",
        "local": "human_g1k_v37_decoy.fasta.gz",
        "decompressed": "human_g1k_v37_decoy.fasta",
        "description": "Reference genome FASTA file",
    },
    "fai": {
        "remote": "human_g1k_v37_decoy.fasta.fai.gz",
        "local": "human_g1k_v37_decoy.fasta.fai.gz",
        "decompressed": "human_g1k_v37_decoy.fasta.fai",
        "description": "FASTA index file",
    },
    "dict": {
        "remote": "human_g1k_v37_decoy.dict.gz",
        "local": "human_g1k_v37_decoy.dict.gz",
        "decompressed": "human_g1k_v37_decoy.dict",
        "description": "Sequence dictionary file",
    },
}


def calculate_md5(filepath, chunk_size=8192):
    """Calculate MD5 checksum of a file."""
    md5_hash = hashlib.md5()
    with open(filepath, "rb") as f:
        while chunk := f.read(chunk_size):
            md5_hash.update(chunk)
    return md5_hash.hexdigest()


def download_file_ftp(remote_path, local_path, description):
    """Download a file using wget from FTP server."""
    ftp_url = f"ftp://{FTP_USERNAME}@{FTP_SERVER}{FTP_BASE_PATH}/{remote_path}"

    logger.info(f"Downloading {description}: {remote_path}")
    logger.info(f"URL: {ftp_url}")

    # Check if file already exists
    if local_path.exists():
        size_mb = local_path.stat().st_size / (1024 * 1024)
        logger.info(f"{local_path} already exists ({size_mb:.1f} MB), skipping download")
        return True

    # Download with wget
    cmd = ["wget", "-c", "--progress=bar:force", "-O", str(local_path), ftp_url]

    logger.info(f"Running: wget -c -O {local_path.name} [ftp_url]")

    process = subprocess.Popen(cmd, stderr=subprocess.PIPE, universal_newlines=True)

    # Monitor progress
    while True:
        output = process.stderr.readline()
        if output == "" and process.poll() is not None:
            break
        if output and "%" in output:
            print(f"\r{local_path.name}: {output.strip()}", end="", flush=True)

    print()  # New line after progress

    if process.returncode != 0:
        logger.error(f"Failed to download {remote_path}")
        return False

    # Report download success
    size_mb = local_path.stat().st_size / (1024 * 1024)
    logger.info(f"Downloaded {local_path.name} successfully ({size_mb:.1f} MB)")

    return True


def decompress_gzip(gz_file, output_file, keep_compressed=True):
    """Decompress a gzip file."""
    logger.info(f"Decompressing {gz_file.name} -> {output_file.name}")

    # Check if decompressed file already exists
    if output_file.exists():
        size_mb = output_file.stat().st_size / (1024 * 1024)
        logger.info(f"{output_file} already exists ({size_mb:.1f} MB), skipping decompression")
        return True

    try:
        with gzip.open(gz_file, "rb") as f_in:
            with open(output_file, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)

        size_mb = output_file.stat().st_size / (1024 * 1024)
        logger.info(f"Decompressed to {output_file.name} ({size_mb:.1f} MB)")

        # Optionally remove compressed file
        if not keep_compressed:
            logger.info(f"Removing compressed file: {gz_file}")
            gz_file.unlink()

        return True

    except Exception as e:
        logger.error(f"Failed to decompress {gz_file}: {e}")
        return False


def verify_reference_integrity(fasta_file, fai_file, dict_file):
    """Perform basic integrity checks for reference files."""
    logger.info("Performing integrity checks...")

    # Check FASTA file exists and is not empty
    if not fasta_file.exists():
        logger.error(f"FASTA file not found: {fasta_file}")
        return False

    fasta_size = fasta_file.stat().st_size
    if fasta_size == 0:
        logger.error(f"FASTA file is empty: {fasta_file}")
        return False

    # Check FAI file exists and is not empty
    if not fai_file.exists():
        logger.error(f"Index file not found: {fai_file}")
        return False

    fai_size = fai_file.stat().st_size
    if fai_size == 0:
        logger.error(f"Index file is empty: {fai_file}")
        return False

    # Check dictionary file exists and is not empty
    if not dict_file.exists():
        logger.error(f"Dictionary file not found: {dict_file}")
        return False

    dict_size = dict_file.stat().st_size
    if dict_size == 0:
        logger.error(f"Dictionary file is empty: {dict_file}")
        return False

    # Basic file size checks (human_g1k_v37_decoy should be ~3GB)
    fasta_size_gb = fasta_size / (1024 * 1024 * 1024)
    if fasta_size_gb < 2.5:  # Should be around 3GB
        logger.warning(f"FASTA file seems small ({fasta_size_gb:.1f} GB), expected ~3GB")

    # Check FAI has reasonable number of sequences
    with open(fai_file) as f:
        fai_lines = f.readlines()

    if len(fai_lines) < 20:  # Should have 25+ chromosomes/contigs including decoys
        logger.warning(f"Index file has only {len(fai_lines)} sequences, expected 25+")

    # Check dictionary has reasonable number of sequences
    with open(dict_file) as f:
        dict_lines = [line for line in f if line.startswith("@SQ")]

    if len(dict_lines) < 20:  # Should have 25+ @SQ lines
        logger.warning(f"Dictionary file has only {len(dict_lines)} sequences, expected 25+")

    # Cross-check that FAI and dictionary have same number of sequences
    if len(fai_lines) != len(dict_lines):
        logger.warning(f"Sequence count mismatch: FAI={len(fai_lines)}, dict={len(dict_lines)}")

    logger.info("Integrity checks passed:")
    logger.info(f"  - FASTA: {fasta_size_gb:.2f} GB")
    logger.info(f"  - Index: {len(fai_lines)} sequences")
    logger.info(f"  - Dictionary: {len(dict_lines)} sequences")

    return True


def main():
    """Download human_g1k_v37_decoy reference genome from Broad Institute FTP."""
    parser = argparse.ArgumentParser(
        description="Download human_g1k_v37_decoy reference genome from Broad Institute FTP"
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("reference"),
        help="Output directory for reference files",
    )
    parser.add_argument(
        "--keep-compressed",
        action="store_true",
        help="Keep compressed .gz files after decompression",
    )
    parser.add_argument(
        "--compressed-only",
        action="store_true",
        help="Download compressed files only (do not decompress)",
    )

    args = parser.parse_args()

    logger.info(f"Starting reference genome download - Log file: {log_file}")
    logger.info(f"Output directory: {args.output_dir}")

    # Check dependencies
    try:
        subprocess.run(["wget", "--version"], capture_output=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        logger.error("wget not found. Please install wget first.")
        sys.exit(1)

    # Create output directory
    args.output_dir.mkdir(parents=True, exist_ok=True)

    # Download files
    logger.info("Downloading reference genome files from Broad Institute FTP...")

    success_count = 0
    for file_type, file_info in REFERENCE_FILES.items():
        remote_name = file_info["remote"]
        local_path = args.output_dir / file_info["local"]
        description = file_info["description"]

        if download_file_ftp(remote_name, local_path, description):
            success_count += 1
        else:
            logger.error(f"Failed to download {file_type} file")

    if success_count != len(REFERENCE_FILES):
        logger.error(f"Only {success_count}/{len(REFERENCE_FILES)} files downloaded successfully")
        sys.exit(1)

    # Decompress files if requested
    if not args.compressed_only:
        logger.info("Decompressing files...")

        for file_type, file_info in REFERENCE_FILES.items():
            gz_path = args.output_dir / file_info["local"]
            output_path = args.output_dir / file_info["decompressed"]

            if not decompress_gzip(gz_path, output_path, args.keep_compressed):
                logger.error(f"Failed to decompress {file_type} file")
                sys.exit(1)

        # Verify integrity of decompressed files
        fasta_file = args.output_dir / REFERENCE_FILES["fasta"]["decompressed"]
        fai_file = args.output_dir / REFERENCE_FILES["fai"]["decompressed"]
        dict_file = args.output_dir / REFERENCE_FILES["dict"]["decompressed"]

        if not verify_reference_integrity(fasta_file, fai_file, dict_file):
            logger.error("Reference integrity checks failed")
            sys.exit(1)

    # Final summary
    logger.info("\n" + "=" * 60)
    logger.info("DOWNLOAD SUMMARY")
    logger.info("=" * 60)

    # List downloaded files
    downloaded_files = list(args.output_dir.glob("*"))
    if downloaded_files:
        logger.info(f"Successfully downloaded {len(downloaded_files)} files:")
        total_size = 0
        for file_path in sorted(downloaded_files):
            if file_path.is_file():
                size_mb = file_path.stat().st_size / (1024 * 1024)
                total_size += size_mb
                logger.info(f"  - {file_path.name}: {size_mb:.1f} MB")
        logger.info(f"Total size: {total_size:.1f} MB")

    logger.info(f"\nReference files available in: {args.output_dir}")

    if not args.compressed_only:
        logger.info("\nTo use with variant calling tools:")
        fasta_name = REFERENCE_FILES["fasta"]["decompressed"]
        logger.info(f"  Reference FASTA: {args.output_dir / fasta_name}")
        logger.info(f"  Index file: {args.output_dir / REFERENCE_FILES['fai']['decompressed']}")
        logger.info(
            f"  Dictionary file: {args.output_dir / REFERENCE_FILES['dict']['decompressed']}"
        )

    logger.info(f"\nLog saved to: {log_file}")
    logger.info("Done!")


if __name__ == "__main__":
    main()
