# VariantCentrifuge

**VariantCentrifuge** is a Python-based command-line tool designed to filter, extract, and refine genetic variant data (VCF files) based on genes of interest, rarity criteria, and impact annotations. Originally inspired by a Bash/R pipeline, VariantCentrifuge provides a more modular, maintainable, and extensible codebase, making it easier to incorporate additional filtering steps, analysis, and reporting formats.

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
│  ├─ config.json
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
  - `pandas` (for XLSX conversion)
  - `pytest` (for testing)
  - `scipy` (for Fisher exact test in variant analysis)
  - Any others specified in `requirements.txt`.

## Installation

1. **Clone the repository:**
   ```sh
   git clone https://github.com/yourusername/VariantCentrifuge.git
   cd VariantCentrifuge
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
   Ensure `bcftools`, `snpEff`, `SnpSift`, and `bedtools` are installed and available in your `PATH` (as described above).

## Usage

**Basic command:**

```sh
variantcentrifuge \
  --gene-name BICC1 \
  --vcf-file path/to/your.vcf \
  --output-file output.tsv
```

**Additional options:**

- `--config CONFIG_FILE` to load default parameters from a JSON config file.
- `--reference REFERENCE` to specify the snpEff reference database.
- `--filters "FILTER_EXPRESSION"` to apply custom SnpSift filters.
- `--fields "FIELD_LIST"` to extract custom fields from the VCF.
- `--samples-file SAMPLES.TXT` to specify the file used in genotype replacement.
- `--xlsx` to convert the final output TSV into XLSX format.
- `--version` to show the current version and exit.

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

The tool loads defaults from a JSON config file (default `config.json`), which can be overridden by `--config`.

## Testing

Run tests with:

```sh
pytest tests/
```

## Roadmap

- **Implement genotype replacement logic** in `replacer.py` more fully.
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
