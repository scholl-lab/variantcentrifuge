# 1000 Genomes two-cohort fixture

Real genotypes for an end-to-end test of association analysis, `--score-artifact` and
`variantcentrifuge meta` (issue #112; first real-data dataset for #60).

| File | Content |
|---|---|
| `1000g_nephro.GRCh38.annotated.vcf.gz` | 666 coding/splice variants in 10 nephrology genes (PKD1, PKD2, PKHD1, COL4A3, COL4A4, COL4A5, HNF1B, UMOD, NPHS1, NPHS2), 1039 samples, GT only, INFO AC/AN/AF (EUR+SAS), snpEff `ANN`/`LOF`/`NMD`, ClinVar `CLNSIG`/`CLNREVSTAT` |
| `populations.tsv` | sample → superpopulation; the two cohorts are EUR (525) and SAS (514) unrelated founders |
| `EUR.cases.txt`, `EUR.controls.txt`, `SAS.cases.txt`, `SAS.controls.txt` | seeded synthetic phenotype |
| `genes.GRCh38.bed` | gene regions from the snpEff GRCh38.p14 database (chr-prefixed) |
| `PROVENANCE.txt`, `SHA256SUMS` | input URLs, tool versions, checksums |

## Sources

- Genotypes: 1000 Genomes 30x on GRCh38, 2022-04-22 phased SNV/INDEL/SV panel (IGSR, open
  access). chrX is the `v2` release; males are haploid outside the PARs.
- Annotation matches the AGDE cohort (`multisample-vcf-cohorts`, `sm-vcf-annotation`):
  snpEff 5.2 `-lof` with the `GRCh38.p14` database, ClinVar weekly release 2025-05-04
  (NCBI archive, public domain). dbNSFP and HGMD are left out because their licences do
  not allow redistribution; nothing in the test needs them.
- Phenotype (`scripts/make_1000g_phenotypes.py`, seed 112): carriers of a rare (AF < 1%)
  HIGH or MODERATE PKHD1 variant are cases with probability 0.6, everyone else 0.2. PKHD1
  is the only true signal; the other genes are real-genotype nulls. The labels are
  synthetic and say nothing about these individuals.

## Rebuild

Needs Docker, the variantcentrifuge image and the snpEff `GRCh38.p14` database
(`snpEff download GRCh38.p14`, 638 MB):

```bash
make fixture-1000g SNPEFF_DATA=/path/to/snpeff/data E2E_IMAGE=variantcentrifuge:0.18.0
cd tests/fixtures/1000g_cohort && sha256sum -c SHA256SUMS
```

The build writes no timestamps or temp paths, so a rebuild with the same image is
byte-identical. Run the test with `make test-e2e` (about 2 minutes).
