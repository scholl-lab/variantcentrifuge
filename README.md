# VariantCentrifuge

**VariantCentrifuge** is a Python-based command-line tool designed to filter, extract, and refine genetic variant data (VCF files) based on genes of interest, rarity criteria, and impact annotations. Built with modularity and extensibility in mind, VariantCentrifuge replaces the complexity of traditional Bash/R pipelines with a cleaner, maintainable Python codebase.

## Key Features

- **Gene-Centric Filtering:**  
  Extract variants from regions defined by genes of interest, using `snpEff` genes2bed to generate BED files.

- **Rare Variant Identification:**  
  Apply custom filters via `SnpSift` to isolate rare and moderate/high-impact variants.

- **Flexible Field Extraction:**  
  Easily specify which fields to extract from the VCF (e.g., gene annotations, functional predictions, allele counts).

- **Genotype Replacement:**  
  Replace genotype fields with corresponding sample IDs, enabling more interpretable variant reports.

- **Phenotype Integration:**  
  Integrate phenotype data from a provided table (CSV or TSV) to further filter or annotate variants based on sample-level attributes.

- **Variant and Gene-Level Analysis:**  
  Perform gene burden analyses (e.g., Fisher’s exact test) and variant-level statistics.  

- **Reporting and Visualization:**  
  - Generate tab-delimited outputs by default and optionally convert them into Excel (XLSX) format.
  - Create an interactive HTML report with sortable variant tables and IGV.js integration for genomic visualization.

## Project Structure

A typical directory layout is:

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
│  ├─ gene_burden.py
│  ├─ generate_html_report.py
│  ├─ generate_igv_report.py
│  ├─ helpers.py
│  ├─ phenotype_filter.py
│  ├─ phenotype.py
│  ├─ pipeline.py
│  ├─ replacer.py
│  ├─ stats.py
│  ├─ utils.py
│  ├─ validators.py
│  └─ templates/
│     └─ index.html
├─ tests/
│  ├─ test_cli.py
│  └─ test_filters.py
├─ requirements.txt
├─ setup.py
├─ pyproject.toml
├─ MANIFEST.in
├─ README.md
└─ LICENSE
```

## Dependencies

- **Python 3.7+**

- **External Tools:**  
  - `snpEff` for generating gene BED files and functional annotations.
  - `SnpSift` for filtering and field extraction.
  - `bcftools` for variant extraction and manipulation.
  - `bedtools` (specifically `sortBed`) for sorting BED files.
  
  **Installation via mamba/conda:**
  ```sh
  mamba create -y -n annotation bcftools snpsift snpeff bedtools
  mamba activate annotation
  ```
  
  Ensure these tools are in your `PATH` before running VariantCentrifuge.

- **Python Packages:**  
  The required Python packages can be installed via `pip` or `mamba/conda`.  
  Minimal required packages include:
  - `pandas` (for XLSX conversion and data handling)
  - `pytest` (for testing)
  - `scipy` (for Fisher exact test in variant analysis)
  - `statsmodels` (for multiple testing correction in gene burden analysis)
  - `jinja2` (for HTML template rendering)
  - `openpyxl` (for XLSX creation)

  To install using `pip`:
  ```sh
  pip install -r requirements.txt
  ```

  Or using `mamba/conda`:
  ```sh
  mamba install pandas pytest scipy statsmodels jinja2 openpyxl
  ```

## Installation

1. **Clone the repository:**
   ```sh
   git clone https://github.com/scholl-lab/variantcentrifuge/.git
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

4. **Check external tools:**
   Ensure `bcftools`, `snpEff`, `SnpSift`, and `bedtools` are installed and available in your `PATH`.

## Configuration

VariantCentrifuge uses a JSON configuration file (`config.json`) to set default parameters. You can specify a custom configuration file with `--config`. If no configuration file is found, a helpful error message will guide you to create one.

**Required Keys:**
- **reference** (`str`): Reference genome database for `snpEff`. No default; must be provided.
- **filters** (`str`): A SnpSift filter expression to select variants. No default; must be provided.
- **fields_to_extract** (`str`): Space-separated list of fields to extract via SnpSift. No default; must be provided.

**Optional Keys and Their Defaults:**
- **interval_expand** (`int`): Number of bases to expand around genes. *Default: 0*  
- **add_chr** (`bool`): Add "chr" prefix to chromosome names. *Default: true*  
- **debug_level** (`str`): Logging level: "DEBUG", "INFO", "WARN", "ERROR". *Default: "INFO"*  
- **no_stats** (`bool`): Skip statistics computation. *Default: false*  
- **perform_gene_burden** (`bool`): Perform gene burden analysis. *Default: false*  
- **gene_burden_mode** (`str`): "samples" or "alleles". *Default: "alleles"*  
- **correction_method** (`str`): "fdr" or "bonferroni" for multiple testing correction. *Default: "fdr"*  
- **igv_enabled** (`bool`): Enable IGV.js integration. *Default: false*  
- **bam_mapping_file** (`str`): Required if igv_enabled=true. No default.  
- **igv_reference** (`str`): Required if igv_enabled=true. No default.

**Example `config.json`:**
```json
{
  "reference": "GRCh37.75",
  "filters": "(( dbNSFP_gnomAD_exomes_AC[0] <= 2 ) | ( na dbNSFP_gnomAD_exomes_AC[0] )) & ((ANN[ANY].IMPACT has 'HIGH') | (ANN[ANY].IMPACT has 'MODERATE'))",
  "fields_to_extract": "CHROM POS REF ALT ANN[0].GENE ANN[0].IMPACT GEN[*].GT",
  "interval_expand": 0,
  "add_chr": true,
  "debug_level": "INFO",
  "no_stats": false,
  "perform_gene_burden": false,
  "gene_burden_mode": "alleles",
  "correction_method": "fdr",
  "igv_enabled": false
}
```

If `config.json` is missing or incomplete, VariantCentrifuge will print a clear error message. Provide required keys in the config or use CLI arguments to override defaults. This encourages a user-friendly configuration workflow.

## Usage

**Basic command:**
```sh
variantcentrifuge \
  --gene-name BICC1 \
  --vcf-file path/to/your.vcf \
  --output-file output.tsv
```

**Additional options:**
- `--config CONFIG_FILE` to load custom parameters from a JSON config file.
- `--reference REFERENCE` to specify the snpEff reference database (overrides config).
- `--filters "FILTER_EXPRESSION"` to apply custom SnpSift filters (overrides config).
- `--fields "FIELD_LIST"` to extract custom fields from the VCF (overrides config).
- `--gene-file GENES.TXT` to provide multiple genes of interest.
- `--samples-file SAMPLES.TXT` for genotype replacement mapping.
- `--phenotype-file PHENO.TSV` along with `--phenotype-sample-column` and `--phenotype-value-column`.
- `--xlsx` to convert the final output TSV into XLSX format.
- `--perform-gene-burden` to run gene burden analysis.
- `--html-report` to generate an interactive HTML report.
- `--igv` with `--bam-mapping-file` and `--igv-reference` for IGV.js integration.
- `--version` to show the current version and exit.

**Example:**
```sh
variantcentrifuge \
  --gene-name BICC1 \
  --vcf-file input.vcf.gz \
  --filters "(( dbNSFP_gnomAD_exomes_AC[0] <= 2 ) | ( na dbNSFP_gnomAD_exomes_AC[0] )) & ((ANN[ANY].IMPACT has 'HIGH') | (ANN[ANY].IMPACT has 'MODERATE'))" \
  --xlsx
```

## Phenotype Integration

If you provide a `--phenotype-file` (CSV or TSV) along with `--phenotype-sample-column` and `--phenotype-value-column`, VariantCentrifuge will integrate sample phenotypes into the final output. This enables downstream filtering or annotation by phenotype.

## Testing

Run tests with:
```sh
pytest tests/
```

## Contributing

Contributions are welcome! Open issues, submit pull requests, or suggest features. Please maintain code quality, follow PEP8 style guidelines, and ensure that all tests pass before submitting a pull request.

## License

This project is licensed under the [MIT License](LICENSE).

## Acknowledgments

- Inspired by prior Bash/R pipelines for variant filtering.
- Built upon the rich ecosystem of bioinformatics tools (snpEff, SnpSift, bcftools, bedtools).
- Special thanks to contributors and the open-source community.
