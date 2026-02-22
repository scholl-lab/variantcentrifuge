# Usage Guide

## Basic Usage

The most basic command to run VariantCentrifuge:

```bash
variantcentrifuge \
  --gene-name BICC1 \
  --vcf-file path/to/your.vcf \
  --output-file output.tsv
```

## Pipeline Flow

VariantCentrifuge processes variants through a multi-stage pipeline. The stages executed depend on the options you provide:

```{mermaid}
flowchart TD
    VCF["VCF File"] --> BED["Gene BED Creation"]
    BED --> PRE{"bcftools prefilter?"}
    PRE -->|yes| BCF["bcftools Prefilter"]
    PRE -->|no| EXT["Variant Extraction"]
    BCF --> EXT
    EXT --> FILT["SnpSift Filter"]
    FILT --> SPLIT{"Split snpEff lines?"}
    SPLIT -->|yes| SPL["Annotation Splitter"]
    SPLIT -->|no| FE["Field Extraction"]
    SPL --> FE
    FE --> SORT["Data Sorting"]
    SORT --> GT{"Genotype replacement?"}
    GT -->|yes| REP["Genotype Replacement"]
    GT -->|no| DF["DataFrame Loading"]
    REP --> DF
    DF --> INH{"Inheritance analysis?"}
    INH -->|yes| PED["Inheritance Analysis"]
    INH -->|no| ANN["Custom Annotations"]
    PED --> ANN
    ANN --> SCR{"Scoring model?"}
    SCR -->|yes| SCORE["Variant Scoring"]
    SCR -->|no| STAT["Statistics"]
    SCORE --> STAT
    STAT --> GB{"Gene burden?"}
    GB -->|yes| BURDEN["Gene Burden Analysis"]
    GB -->|no| FF{"Final filter?"}
    BURDEN --> FF
    FF -->|yes| FFILT["Final Filter"]
    FF -->|no| PSE{"Pseudonymize?"}
    FFILT --> PSE
    PSE -->|yes| PSEUDO["Pseudonymization"]
    PSE -->|no| OUT["TSV Output"]
    PSEUDO --> OUT
    OUT --> XLSX{"Excel?"}
    OUT --> HTML{"HTML report?"}
    OUT --> IGV{"IGV report?"}
    XLSX -->|yes| XL["Excel Report"]
    HTML -->|yes| HR["HTML Report"]
    IGV -->|yes| IR["IGV Report"]

    style VCF fill:#4a90d9,color:#fff
    style OUT fill:#2ecc71,color:#fff
    style XL fill:#2ecc71,color:#fff
    style HR fill:#2ecc71,color:#fff
    style IR fill:#2ecc71,color:#fff
```

## Command Line Options

### General Options

| Flag | Default | Description |
|------|---------|-------------|
| `--version` | — | Show version and exit |
| `--log-level LEVEL` | `INFO` | Logging level: `DEBUG`, `INFO`, `WARN`, `ERROR` |
| `--log-file PATH` | — | Write logs to file in addition to stderr |
| `-c` / `--config PATH` | — | Path to JSON configuration file |

### Core Input/Output

| Flag | Default | Description |
|------|---------|-------------|
| `-v` / `--vcf-file PATH` | **required** | Input VCF file (can be gzip-compressed) |
| `-o` / `--output-file PATH` | — | Output TSV file path, or `stdout`/`-` for stdout |
| `--output-dir DIR` | `output` | Directory for intermediate and final output |
| `--xlsx` | `false` | Also produce Excel output |
| `--keep-intermediates` | `false` | Keep all intermediate files after run |
| `--archive-results` | `false` | Create timestamped `.tar.gz` archive of results |

### Gene Selection

| Flag | Default | Description |
|------|---------|-------------|
| `-g` / `--gene-name GENES` | — | Gene name(s) (space-separated string) |
| `-G` / `--gene-file PATH` | — | File with gene names, one per line |
| `--transcript-list IDS` | — | Comma-separated transcript IDs to filter for (e.g., `NM_007294.4,NM_000059.4`) |
| `--transcript-file PATH` | — | File with transcript IDs, one per line |

### Filtering Options

| Flag | Default | Description |
|------|---------|-------------|
| `-r` / `--reference DB` | from config | snpEff reference database (e.g., `GRCh37.p13`) |
| `-f` / `--filters EXPR` | from config | SnpSift filter expression |
| `--bcftools-prefilter EXPR` | — | bcftools expression for early variant pre-filtering during extraction |
| `--preset NAME` | — | Apply named filter preset from config (repeatable, combined with AND). Use commas within a single `--preset` for multiple: `--preset rare,coding` |
| `--field-profile PROFILE` | `dbnsfp4` | Field profile for annotation database compatibility (e.g., `dbnsfp4`, `dbnsfp5`). See [Field Profiles](configuration.md#field-profiles) |
| `--list-field-profiles` | — | List available field profiles and exit |
| `--late-filtering` | `false` | Apply SnpSift filters after scoring (allows filtering on computed scores) |
| `--final-filter EXPR` | — | Pandas `query()` expression applied after all annotations and scores |
| `--split-snpeff-lines MODE` | — | Split multi-annotation lines: `before_filters` or `after_filters`. Omit to skip |

### Tumor-Normal Filtering

These flags configure somatic variant presets (`somatic`, `loh`, `tumor_only`) by expanding template variables:

| Flag | Default | Description |
|------|---------|-------------|
| `--tumor-sample-index N` | `1` | 0-based index of tumor sample in VCF |
| `--normal-sample-index N` | `0` | 0-based index of normal sample in VCF |
| `--tumor-dp-min N` | `20` | Minimum read depth for tumor sample |
| `--normal-dp-min N` | `20` | Minimum read depth for normal sample |
| `--tumor-af-min F` | `0.05` | Minimum allele frequency for tumor sample |
| `--normal-af-max F` | `0.03` | Maximum allele frequency for normal sample |

### VCF Annotation Inspection

| Flag | Default | Description |
|------|---------|-------------|
| `--show-vcf-annotations` | `false` | Print VCF INFO/FORMAT fields and exit |
| `--annotation-filter TEXT` | — | Filter annotation output by substring (case-insensitive) |
| `--annotation-format FMT` | `table` | Output format: `table` or `json` |
| `--info-only` | `false` | Show only INFO fields (mutually exclusive with `--format-only`) |
| `--format-only` | `false` | Show only FORMAT fields (mutually exclusive with `--info-only`) |

### Field Extraction & Formatting

| Flag | Default | Description |
|------|---------|-------------|
| `-e` / `--fields FIELDS` | from config | Fields to extract with SnpSift extractFields |
| `--no-replacement` | `false` | Skip genotype replacement step |
| `--add-column NAME` | — | Add blank column(s) to final output (repeatable) |
| `--no-links` | `false` | Disable adding URL link columns |
| `--add-chr` | `false` | Add `chr` prefix to chromosome names in BED |
| `--remove-sample-substring TEXT` | — | Strip this substring from all sample names |

### Genotype Analysis

| Flag | Default | Description |
|------|---------|-------------|
| `--genotype-filter FILTER` | — | Genotype filter: `het`, `hom`, `comp_het`, or comma-separated (e.g., `het,comp_het`) |
| `--gene-genotype-file PATH` | — | TSV with per-gene genotype rules (columns: `GENE`, `GENOTYPES`). Overrides global `--genotype-filter` |
| `--append-extra-sample-fields [FIELDS...]` | — | Append extra FORMAT fields (e.g., `DP AD`) next to genotypes |
| `--extra-sample-field-delimiter CHAR` | `:` | Delimiter between genotype and extra fields |

### Inheritance Analysis

| Flag | Default | Description |
|------|---------|-------------|
| `--ped PATH` | — | PED file defining family structure for inheritance analysis |
| `--inheritance-mode MODE` | `simple` | Output format: `simple` (pattern only), `columns` (separate columns), `full` (complete JSON) |
| `--no-vectorized-comp-het` | `false` | Use original compound het implementation instead of vectorized (10-50x slower) |

Inheritance analysis is triggered when `--ped` or `--inheritance-mode` is specified. Without `--ped`, all samples are treated as affected singletons.

**Supported inheritance patterns:** de novo, autosomal dominant (AD), autosomal recessive (AR), X-linked recessive (XLR), X-linked dominant (XLD), compound heterozygous, mitochondrial.

### ClinVar PM5 Annotation

| Flag | Default | Description |
|------|---------|-------------|
| `--clinvar-pm5-lookup PATH` | — | Path to pre-built PM5 lookup table. Adds `PM5`, `PM5_evidence_count`, and `PM5_known_variants` columns |

### Phenotype & Sample Groups

| Flag | Default | Description |
|------|---------|-------------|
| `--phenotype-file PATH` | — | Path to phenotype file (`.csv` or `.tsv`) |
| `--phenotype-sample-column NAME` | — | Column name for sample IDs in phenotype file |
| `--phenotype-value-column NAME` | — | Column name for phenotype values |
| `--case-phenotypes TERMS` | — | Comma-separated HPO terms defining case group |
| `--control-phenotypes TERMS` | — | Comma-separated HPO terms defining control group |
| `--case-phenotypes-file PATH` | — | File with HPO terms for case group |
| `--control-phenotypes-file PATH` | — | File with HPO terms for control group |
| `--case-samples IDS` | — | Comma-separated sample IDs for case group |
| `--control-samples IDS` | — | Comma-separated sample IDs for control group |
| `--case-samples-file PATH` | — | File with sample IDs for case group |
| `--control-samples-file PATH` | — | File with sample IDs for control group |

### Statistical Analysis

| Flag | Default | Description |
|------|---------|-------------|
| `--perform-gene-burden` | `false` | Run gene burden analysis |
| `--gene-burden-mode MODE` | `alleles` | Gene burden mode: `samples` or `alleles` |
| `--correction-method METHOD` | `fdr` | Multiple testing correction: `fdr` or `bonferroni` |
| `--no-stats` | `false` | Skip statistics computation step |
| `--stats-output-file PATH` | — | File to write analysis statistics |
| `--stats-config PATH` | — | Path to custom statistics configuration JSON |

### Association Testing

For modular rare variant association testing beyond the basic Fisher's exact test in `--perform-gene-burden`:

| Flag | Default | Description |
|------|---------|-------------|
| `--perform-association` | `false` | Run the modular association testing framework |
| `--association-tests TESTS` | `fisher` | Comma-separated tests: `fisher`, `logistic_burden`, `linear_burden`, `skat`, `skat_o`, `coast` |
| `--skat-backend BACKEND` | `python` | SKAT backend: `python` (default, thread-safe) or `r` (deprecated) |
| `--coast-backend BACKEND` | `python` | COAST backend: `python` (default) or `r` (deprecated) |
| `--covariate-file PATH` | — | TSV/CSV covariate file (first column = sample ID, header required) |
| `--covariates NAMES` | — | Comma-separated covariate column names to include (default: all columns) |
| `--categorical-covariates NAMES` | — | Comma-separated columns to treat as categorical (auto-detected otherwise) |
| `--trait-type TYPE` | `binary` | Trait type for burden tests: `binary` or `quantitative` |
| `--pca-file PATH` | — | PCA file (PLINK `.eigenvec`, AKT output, or generic TSV) |
| `--pca-tool TOOL` | — | Set to `akt` to invoke AKT subprocess for PCA computation |
| `--pca-components N` | `10` | Number of PCA components to include as covariates |
| `--variant-weights SCHEME` | `beta:1,25` | Variant weight scheme: `beta:a,b`, `uniform`, `cadd`, `revel`, `combined` |
| `--variant-weight-params JSON` | — | JSON string of weight scheme parameters (e.g., `'{"alpha": 1, "beta": 25}'`) |
| `--coast-weights WEIGHTS` | — | COAST category weights as comma-separated floats (BMV,DMV,PTV) |
| `--diagnostics-output DIR` | — | Directory for diagnostics output: `lambda_gc.tsv`, `qq_data.tsv`, `summary.txt` |

For detailed usage, test selection guidance, and examples, see the [Association Testing Guide](guides/association_testing.md).

### Scoring & Custom Annotations

| Flag | Default | Description |
|------|---------|-------------|
| `--scoring-config-path DIR` | — | Directory containing scoring model (`variable_assignment_config.json` + `formula_config.json`) |
| `--annotate-bed PATH` | — | BED file for region annotation (repeatable) |
| `--annotate-gene-list PATH` | — | Gene list file — adds yes/no column per file (repeatable) |
| `--annotate-json-genes PATH` | — | JSON gene data file (repeatable, requires `--json-gene-mapping`) |
| `--json-gene-mapping JSON` | — | JSON string mapping fields: `'{"identifier":"gene_symbol","dataFields":["panel","inheritance"]}'` |
| `--json-genes-as-columns` | `false` | Add each dataField as its own column instead of `Custom_Annotation` |

### Reporting & Visualization

| Flag | Default | Description |
|------|---------|-------------|
| `--html-report` | `false` | Generate interactive HTML report with filtering, charts, and summary dashboard |
| `--igv` | `false` | Enable IGV.js integration (requires `--bam-mapping-file` and a genome reference) |
| `--bam-mapping-file PATH` | — | TSV/CSV mapping sample IDs to BAM files |
| `--igv-reference REF` | — | Genome reference for IGV (e.g., `hg19`, `hg38`) |
| `--igv-fasta PATH` | — | Local FASTA file for IGV (overrides `--igv-reference`) |
| `--igv-ideogram PATH` | — | Ideogram file for chromosome visualization |
| `--igv-flanking N` | `50` | Flanking region in bp for IGV reports |

### Performance & Processing

| Flag | Default | Description |
|------|---------|-------------|
| `--threads N` | `auto` | CPU cores for parallel processing (`auto` detects available cores) |
| `--no-chunked-processing` | `false` | Disable chunked processing (may cause memory issues on large files) |
| `--force-chunked-processing` | `false` | Force chunked processing even for small files |
| `--sort-memory-limit SIZE` | `auto` | Memory for external sort (e.g., `4G`, `auto` = 10% of available) |
| `--sort-parallel N` | `4` | Parallel threads for sorting |
| `--genotype-replacement-method METHOD` | `auto` | Method: `auto`, `sequential`, `vectorized`, `chunked-vectorized`, `parallel`, `streaming-parallel` |
| `--max-memory-gb N` | auto-detected | Maximum memory in GB for processing |
| `--force-inheritance-processing` | `false` | Force inheritance analysis even if exceeding memory limits |
| `--sample-column-creation-method METHOD` | `auto` | Column creation: `iterative`, `vectorized`, or `auto` |
| `--memory-safety-factor F` | `0.92` | Fraction of allocated memory to use (0–1) |
| `--inheritance-memory-fraction F` | `0.85` | Fraction of safe memory for inheritance analysis |

### Checkpoint & Resume

| Flag | Default | Description |
|------|---------|-------------|
| `--enable-checkpoint` | `false` | Enable checkpoint tracking for pipeline state |
| `--resume` | `false` | Resume from last successful checkpoint |
| `--resume-from STAGE` | — | Restart from a specific pipeline stage |
| `--checkpoint-checksum` | `false` | Use SHA256 checksums for file validation (recommended for production) |
| `--show-checkpoint-status` | — | Show checkpoint status and exit |
| `--list-stages` | — | List all available stages for current configuration and exit |
| `--list-checkpoints` | — | List completed stages from checkpoint file and exit |
| `--interactive-resume` | — | Interactively select resume point |

See the [Resume System](resume_system.md) documentation for details.

### Data Privacy

| Flag | Default | Description |
|------|---------|-------------|
| `--pseudonymize` | `false` | Enable sample pseudonymization |
| `--pseudonymize-schema SCHEMA` | `sequential` | Schema: `sequential`, `categorical`, `anonymous`, `custom` |
| `--pseudonymize-prefix TEXT` | `SAMPLE` | Prefix for sequential schema |
| `--pseudonymize-pattern PATTERN` | — | Custom pattern for `custom` schema |
| `--pseudonymize-category-field FIELD` | `phenotype` | Metadata field for categorical schema |
| `--pseudonymize-table PATH` | — | Path to save mapping table (required with `--pseudonymize`) |
| `--pseudonymize-ped` | `false` | Also create pseudonymized PED file |

See [Privacy and Pseudonymization](guides/../user-guide/privacy_and_pseudonymization.md) for details.

### Miscellaneous

| Flag | Default | Description |
|------|---------|-------------|
| `--gzip-intermediates` | `true` | Compress intermediate TSV files (fast level-1 gzip) |
| `--no-gzip-intermediates` | — | Disable intermediate compression |

## Examples

### Basic Gene Analysis

```bash
variantcentrifuge \
  --gene-name BRCA1 \
  --vcf-file samples.vcf.gz \
  --output-file brca1_variants.tsv
```

### Filtered Analysis with Presets

```bash
variantcentrifuge \
  --gene-file cancer_genes.txt \
  --vcf-file samples.vcf.gz \
  --preset rare,coding \
  --html-report \
  --xlsx \
  --output-file cancer_variants.tsv
```

### Family Trio with Inheritance Analysis

```bash
variantcentrifuge \
  --gene-file disease_genes.txt \
  --vcf-file trio.vcf.gz \
  --ped family.ped \
  --inheritance-mode columns \
  --preset rare,coding \
  --html-report \
  --output-file trio_analysis.tsv
```

### Comprehensive Analysis with Reports

```bash
variantcentrifuge \
  --gene-name BRCA1 \
  --vcf-file samples.vcf.gz \
  --phenotype-file patient_data.tsv \
  --phenotype-sample-column "sample_id" \
  --phenotype-value-column "disease_status" \
  --perform-gene-burden \
  --html-report \
  --xlsx \
  --output-file brca1_analysis.tsv
```

### IGV Integration

```bash
variantcentrifuge \
  --gene-name TP53 \
  --vcf-file samples.vcf.gz \
  --igv \
  --bam-mapping-file bam_files.tsv \
  --igv-reference hg38 \
  --html-report \
  --output-file tp53_variants.tsv
```

### Variant Scoring

```bash
variantcentrifuge \
  --gene-file kidney_genes.txt \
  --vcf-file patient.vcf.gz \
  --scoring-config-path scoring/nephro_candidate_score \
  --preset rare,coding \
  --html-report \
  --output-file scored_variants.tsv
```

### Tumor-Normal Somatic Analysis

```bash
variantcentrifuge \
  --gene-file oncogenes.txt \
  --vcf-file tumor_normal.vcf.gz \
  --preset somatic,coding \
  --tumor-sample-index 1 \
  --normal-sample-index 0 \
  --tumor-dp-min 30 \
  --tumor-af-min 0.05 \
  --normal-af-max 0.02 \
  --html-report \
  --output-file somatic_variants.tsv
```

### Custom Annotations

```bash
# Annotate with JSON gene information
variantcentrifuge \
  --gene-name BRCA1 \
  --vcf-file samples.vcf.gz \
  --annotate-json-genes gene_metadata.json \
  --json-gene-mapping '{"identifier":"gene_symbol","dataFields":["panel","inheritance","function"]}' \
  --output-file annotated_variants.tsv

# Multiple annotation sources
variantcentrifuge \
  --gene-file cancer_genes.txt \
  --vcf-file samples.vcf.gz \
  --annotate-bed cancer_hotspots.bed \
  --annotate-gene-list actionable_genes.txt \
  --annotate-json-genes gene_panels.json \
  --json-gene-mapping '{"identifier":"symbol","dataFields":["panel_name","evidence_level"]}' \
  --html-report \
  --output-file multi_annotated.tsv
```

### Advanced Filtering

```bash
# Pre-filter with bcftools for performance
variantcentrifuge \
  --gene-file large_gene_list.txt \
  --vcf-file large_cohort.vcf.gz \
  --bcftools-prefilter 'FILTER="PASS" && INFO/AC<10' \
  --preset rare,coding \
  --output-file filtered_variants.tsv

# Final filter using pandas query syntax
variantcentrifuge \
  --gene-file cancer_genes.txt \
  --vcf-file samples.vcf.gz \
  --preset rare,coding \
  --scoring-config-path scoring/nephro_candidate_score \
  --final-filter 'nephro_candidate_score > 5 and IMPACT == "HIGH"' \
  --output-file high_priority_variants.tsv

# Filter on inheritance patterns
variantcentrifuge \
  --gene-file disease_genes.txt \
  --vcf-file trio.vcf.gz \
  --ped family.ped \
  --inheritance-mode columns \
  --final-filter 'Inheritance_Pattern in ["de_novo", "compound_heterozygous"]' \
  --output-file denovo_and_compound_het.tsv
```

### Checkpoint and Resume

```bash
# Run large analysis with checkpoint tracking
variantcentrifuge \
  --gene-file all_protein_coding_genes.txt \
  --vcf-file large_cohort.vcf.gz \
  --enable-checkpoint \
  --threads 16 \
  --preset rare,coding \
  --html-report \
  --output-file all_genes_analysis.tsv

# If interrupted, resume from last checkpoint
variantcentrifuge \
  --gene-file all_protein_coding_genes.txt \
  --vcf-file large_cohort.vcf.gz \
  --enable-checkpoint \
  --resume \
  --threads 16 \
  --preset rare,coding \
  --html-report \
  --output-file all_genes_analysis.tsv

# Check status of a previous run
variantcentrifuge \
  --show-checkpoint-status \
  --output-dir previous_analysis/
```

### Docker

All examples above work inside the Docker container by mounting your data directory:

```bash
docker run --rm -v ./data:/data \
  ghcr.io/scholl-lab/variantcentrifuge:latest \
  --gene-name BRCA1 \
  --vcf-file /data/input.vcf.gz \
  --preset rare,coding \
  --html-report \
  --output-file /data/output.tsv
```

See the [Installation Guide](installation.md#method-4-docker-recommended-for-quick-setup) for setup details.

## Input File Formats

### VCF Files

- Standard VCF format (v4.0 or later)
- Can be compressed with gzip (.vcf.gz)
- Should be annotated with snpEff for optimal functionality

### Gene Files

Text file with one gene name per line:

```
BRCA1
BRCA2
TP53
ATM
```

### PED Files

Standard PLINK pedigree format (tab-separated, 6 columns):

```
#Family  Individual  Father  Mother  Sex  Affected
FAM001   proband     father  mother  1    2
FAM001   father      0       0       1    1
FAM001   mother      0       0       2    1
```

- Sex: 1=male, 2=female, 0=unknown
- Affected: 1=unaffected, 2=affected, 0=unknown

### Sample Mapping Files

Tab-separated file for genotype replacement:

```
original_id	new_id
sample_001	Patient_A
sample_002	Patient_B
sample_003	Control_001
```

### Phenotype Files

Tab or comma-separated file with sample information:

```
sample_id	disease_status	age	sex
Patient_A	case	45	F
Patient_B	case	52	M
Control_001	control	48	F
```

### BAM Mapping Files

For IGV integration, provide a mapping from sample IDs to BAM file paths:

```
sample_id	bam_path
Patient_A	/path/to/patient_a.bam
Patient_B	/path/to/patient_b.bam
Control_001	/path/to/control_001.bam
```

### JSON Gene Files

For gene annotation, provide a JSON file containing an array of gene objects:

```json
[
  {
    "gene_symbol": "BRCA1",
    "panel": "HereditaryCancer",
    "inheritance": "AD",
    "function": "DNA repair"
  },
  {
    "gene_symbol": "TP53",
    "panel": "HereditaryCancer",
    "inheritance": "AD",
    "function": "Tumor suppressor"
  }
]
```

The `--json-gene-mapping` parameter specifies:
- `identifier`: The field containing the gene symbol (e.g., "gene_symbol")
- `dataFields`: Array of fields to include as annotations (e.g., ["panel", "inheritance", "function"])

## Output Files

### Main Output

- **TSV file** — Tab-separated variant table with extracted fields
- **XLSX file** — Excel format (if `--xlsx` specified)
- **Metadata file** — Analysis parameters and tool versions

### Optional Outputs

- **HTML report** — Interactive variant browser with filtering, summary dashboard, and charts (if `--html-report` specified)
- **IGV reports** — Individual variant visualization (if `--igv` specified)
- **Gene burden results** — Statistical analysis (if `--perform-gene-burden` specified)
- **Pseudonymization mapping** — Sample ID mapping table (if `--pseudonymize` specified)

## Configuration

See the [Configuration Guide](configuration.md) for detailed information about setting up configuration files and customizing VariantCentrifuge behavior.

## Troubleshooting

### Common Issues

1. **No variants found:**
   - Check that your VCF file contains variants in the specified gene regions
   - Verify gene names are correct and match your reference annotation
   - Review filter expressions — they may be too restrictive

2. **External tool errors:**
   - Ensure all required tools are installed and in PATH
   - Check that snpEff database matches your VCF reference
   - Verify file permissions and disk space

3. **Memory issues:**
   - Use `--bcftools-prefilter` to reduce data early
   - Try `--genotype-replacement-method chunked-vectorized`
   - Set `--max-memory-gb` to limit memory usage
   - Use `--threads 1` to reduce concurrent memory pressure

### Getting Help

- Use `variantcentrifuge --help` for command-line options
- Check the [API Reference](api/index.md) for detailed function documentation
- Report issues on [GitHub](https://github.com/scholl-lab/variantcentrifuge/issues)
