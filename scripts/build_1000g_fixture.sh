#!/usr/bin/env bash
# Build tests/fixtures/1000g_cohort: real 1000 Genomes genotypes for 10 nephrology
# genes, annotated like the AGDE cohort (snpEff 5.2 GRCh38.p14 -lof, ClinVar
# 2025-05-04), plus a seeded synthetic phenotype. Run via `make fixture-1000g`,
# which executes this inside the variantcentrifuge image (bcftools, snpEff, SnpSift).
# dbNSFP and HGMD are not added: their licences forbid redistribution.
set -euo pipefail

OUT=${OUT:-tests/fixtures/1000g_cohort}
GENES="PKD1 PKD2 PKHD1 COL4A3 COL4A4 COL4A5 HNF1B UMOD NPHS1 NPHS2"
SNPEFF_DB=GRCh38.p14
KG=https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage
PANEL=$KG/working/20220422_3202_phased_SNV_INDEL_SV/1kGP_high_coverage_Illumina.CHROM.filtered.SNV_INDEL_SV_phased_panel.vcf.gz
PED=$KG/20130606_g1k_3202_samples_ped_population.txt
CLINVAR=https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/2025/clinvar_20250504.vcf.gz

fetch() { python3 -c 'import sys, urllib.request; urllib.request.urlretrieve(*sys.argv[1:])' "$1" "$2"; }

tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT
mkdir -p "$OUT"

# Unrelated founders (no parents in the 3202 pedigree) of two superpopulations = two cohorts.
fetch "$PED" "$tmp/ped.txt"
awk 'NR > 1 && $3 == "0" && $4 == "0" && ($7 == "EUR" || $7 == "SAS") {print $2 "\t" $7}' \
    "$tmp/ped.txt" | sort > "$OUT/populations.tsv"
cut -f1 "$OUT/populations.tsv" > "$tmp/samples.txt"

# Gene regions from the same snpEff database that annotates the variants.
snpEff -Xmx8g genes2bed "$SNPEFF_DB" $GENES 2>/dev/null \
    | awk 'NR > 1 {split($4, a, ";"); print "chr" $1 "\t" $2 "\t" $3 "\t" a[1]}' \
    | sort -k1,1V -k2,2n > "$OUT/genes.GRCh38.bed"

parts=()
for chrom in $(cut -f1 "$OUT/genes.GRCh38.bed" | uniq); do
    grep -P "^$chrom\t" "$OUT/genes.GRCh38.bed" > "$tmp/$chrom.bed"
    url=${PANEL/CHROM/$chrom}
    [ "$chrom" = chrX ] && url=${url/.vcf.gz/.v2.vcf.gz}  # chrX was re-released as v2
    # Run in $tmp: htslib caches the remote .tbi in the working directory.
    (cd "$tmp" && bcftools view --no-version -R "$chrom.bed" -S samples.txt --min-ac 1 -Ou "$url" \
        | bcftools annotate --no-version -x INFO,^FORMAT/GT -Ob -o "$chrom.bcf")
    parts+=("$tmp/$chrom.bcf")
done
bcftools concat --no-version "${parts[@]}" -Ou \
    | bcftools norm --no-version -m -any -Ou \
    | bcftools view --no-version --min-ac 1 -Ou \
    | bcftools +fill-tags --no-version -Ov -o "$tmp/geno.vcf" -- -t AC,AN,AF

# ClinVar subset, renamed to the panel's chr-prefixed contigs.
sed 's/^chr//' "$OUT/genes.GRCh38.bed" > "$tmp/regions_nochr.bed"
for c in $(seq 1 22) X Y MT; do echo -e "$c\tchr$c"; done > "$tmp/chrmap.txt"
fetch "$CLINVAR" "$tmp/clinvar.vcf.gz"
fetch "$CLINVAR.tbi" "$tmp/clinvar.vcf.gz.tbi"
bcftools view --no-version -R "$tmp/regions_nochr.bed" "$tmp/clinvar.vcf.gz" -Ou \
    | bcftools annotate --no-version --rename-chrs "$tmp/chrmap.txt" -Oz -o "$tmp/clinvar.sub.vcf.gz"
tabix -p vcf "$tmp/clinvar.sub.vcf.gz"

# Annotate like AGDE, then drop variants with only MODIFIER consequences (intronic,
# intergenic) to keep the fixture small.
out_abs=$(realpath "$OUT")
cd "$tmp"  # relative paths keep the random temp dir out of the ##SnpEffCmd/##SnpSiftCmd headers
snpEff -Xmx8g ann -noStats -lof "$SNPEFF_DB" geno.vcf \
    | SnpSift -Xmx4g annotate -info CLNSIG,CLNREVSTAT clinvar.sub.vcf.gz - \
    | bcftools view --no-version -i 'INFO/ANN ~ "|HIGH|" || INFO/ANN ~ "|MODERATE|" || INFO/ANN ~ "|LOW|"' \
        -Oz -o "$out_abs/1000g_nephro.GRCh38.annotated.vcf.gz"
cd - > /dev/null
tabix -f -p vcf "$OUT/1000g_nephro.GRCh38.annotated.vcf.gz"

python3 "$(dirname "$0")/make_1000g_phenotypes.py" "$OUT"

{
    echo "panel: ${PANEL/CHROM/<chrom>}"
    echo "pedigree: $PED"
    echo "clinvar: $CLINVAR"
    echo "snpeff_db: $SNPEFF_DB"
    snpEff -version 2>&1 | head -1
    bcftools --version | head -1
} > "$OUT/PROVENANCE.txt"
(cd "$OUT" && sha256sum -- *.vcf.gz *.tsv *.bed *.txt | grep -v SHA256SUMS > SHA256SUMS)
