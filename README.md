# variantcentrifuge

**variantcentrifuge** is a Python-based command-line tool designed to filter, extract, and refine genetic variant data (VCF files) based on genes of interest, rarity criteria, and impact annotations. Originally inspired by a Bash/R pipeline, variantcentrifuge provides a more modular, maintainable, and extensible codebase, making it easier to incorporate additional filtering steps, analysis, and reporting formats.

## Key Features

- **Gene-Centric Filtering:** Extract variants from regions defined by genes of interest, using `snpEff` genes2bed to generate BED files.
- **Rare Variant Identification:** Apply custom filters via `SnpSift` to isolate rare and moderate/high-impact variants.
- **Flexible Field Extraction:** Easily specify which fields to extract from the VCF (e.g., gene annotations, functional predictions, allele counts).
- **Genotype Replacement:** Replace genotype fields with corresponding sample IDs, enabling more interpretable variant reports.
- **Phenotype Integration (Planned):** Integrate phenotype data to further filter variants based on sample-level attributes.
- **Variant Analysis (Planned):** Add optional gene burden or variant-level analyses to derive statistical and clinical insights.
- **Output Formats:** Generate tab-delimited outputs by default and optionally convert them into Excel (XLSX) format for easy downstream use.

## Project Structure

A typical directory layout might look like this:

```
variantcentrifuge/
├─ variantcentrifuge/
│  ├─ __init__.py
│  ├─ analyze_variants.py
│  ├─ cli.py
│  ├─ config.py
│  ├─ converter.py
│  ├─ extractor.py
│  ├─ filters.py
│  ├─ gene_bed.py
│  ├─ phenotype_filter.py
│  ├─ replacer.py
│  └─ utils.py
├─ tests/
│  ├─ test_cli.py
│  └─ test_filters.py
├─ requirements.txt
├─ setup.py
├─ pyproject.toml
├─ README.md
└─ LICENSE
```

## Architecture and Modules

Each step of the pipeline is implemented in a dedicated Python module:

- **`cli.py`**: Handles command-line argument parsing and orchestrates the entire pipeline.
- **`config.py`**: Loads and merges configuration settings from a file and command-line arguments.
- **`gene_bed.py`**: Uses `snpEff` to produce BED files for gene regions, optionally expanding intervals and adding `chr` prefixes.
- **`filters.py`**: Wraps `bcftools` and `SnpSift filter` for extracting and filtering VCF variants.
- **`extractor.py`**: Uses `SnpSift extractFields` to extract the specified fields of interest.
- **`replacer.py`**: Placeholder module for replacing genotype fields with sample IDs.
- **`phenotype_filter.py`**: Placeholder for integrating phenotype-based filtering.
- **`converter.py`**: Converts TSV output to XLSX using `pandas`.
- **`analyze_variants.py`**: Placeholder for gene burden or variant-level analyses.
- **`utils.py`**: Provides utility functions for logging and running external commands.

## Dependencies

- **Python 3.7+**  
- **External Tools:**  
  - `snpEff` for generating gene BED files and functional annotations.
  - `SnpSift` for filtering and field extraction.
  - `bcftools` for variant extraction and manipulation.
  - `bedtools` (specifically `sortBed`) for sorting BED files.
- **Python Packages:**  
  - `pandas` (for XLSX conversion)
  - `pytest` (for testing)
  - Any others specified in `requirements.txt`.

## Installation

1. **Clone the repository:**
   ```sh
   git clone https://github.com/yourusername/variantcentrifuge.git
   cd variantcentrifuge
   ```

2. **Set up a virtual environment (optional but recommended):**
   ```sh
   python3 -m venv venv
   source venv/bin/activate
   ```

3. **Install the tool with pip:**
   ```sh
   pip install .
   ```

   This command will install all Python dependencies listed in `setup.py` and `requirements.txt`, and register a `variantcentrifuge` CLI command.

4. **Check external tools:**
   Ensure `bcftools`, `snpEff`, `SnpSift`, and `bedtools` are installed and available in your `PATH`.

## Usage

**Basic command:**

```sh
variantcentrifuge \
  --gene-name BICC1 \
  --vcf-file path/to/your.vcf \
  --output-file output.tsv
```

Or, if you prefer running via Python directly:
```sh
python -m variantcentrifuge.cli \
  --gene-name BICC1 \
  --vcf-file path/to/your.vcf \
  --output-file output.tsv
```

**Additional options:**

- `--config CONFIG_FILE` to load default parameters from a config file.
- `--reference REFERENCE` to specify the snpEff reference database.
- `--filters "FILTER_EXPRESSION"` to apply custom SnpSift filters.
- `--fields "FIELD_LIST"` to extract custom fields from the VCF.
- `--samples-file SAMPLES.TXT` to specify the file used in genotype replacement.
- `--xlsx` to convert the final output TSV into XLSX format.

## Example

For example, to run a full pipeline extracting variants for gene `BICC1`, filtering for rare and moderate/high-impact variants using default fields, and converting output to XLSX:

```sh
variantcentrifuge \
  --gene-name BICC1 \
  --vcf-file input.vcf.gz \
  --filters "(( dbNSFP_gnomAD_exomes_AC[0] <= 2 ) | ( na dbNSFP_gnomAD_exomes_AC[0] )) & ((ANN[ANY].IMPACT has 'HIGH') | (ANN[ANY].IMPACT has 'MODERATE'))" \
  --xlsx
```

## Configuration

The tool can load defaults from a config file (e.g., `config.cfg`). A typical configuration file might look like:

```ini
[DEFAULT]
reference=GRCh38.mane.1.0.refseq
add_chr=true
filters=(( dbNSFP_gnomAD_exomes_AC[0] <= 2 ) | ( na dbNSFP_gnomAD_exomes_AC[0] )) & ((ANN[ANY].IMPACT has 'HIGH') | (ANN[ANY].IMPACT has 'MODERATE'))
fields_to_extract=CHROM POS REF ALT ID QUAL ...
sample_file=samples.txt
output_file=variants.tsv
```

## Testing

Run tests with:

```sh
pytest tests/
```

Ensure that all external tools are either mocked or installed when running tests.

## Roadmap

- **Implement genotype replacement logic** in `replacer.py`.
- **Integrate phenotype data** for advanced variant filtering.
- **Implement variant/gene burden analysis** in `analyze_variants.py`.
- Add more comprehensive tests and CI/CD support.

## Contributing

Contributions are welcome! Feel free to open issues, submit pull requests, or suggest features. Please maintain code quality, follow PEP8 styling, and ensure that all tests pass.

## License

This project is licensed under the [MIT License](LICENSE).

## Acknowledgments

- Inspired by prior Bash/R pipelines for variant filtering.
- Built upon the rich ecosystem of bioinformatics tools (snpEff, SnpSift, bcftools, bedtools).
- Special thanks to contributors and the open-source community.