#!/usr/bin/env python
"""Create a test VCF file with variants in genes needed for regression testing."""

import gzip
from datetime import datetime

# Test genes needed for regression tests
TEST_GENES = {
    "BRCA1": {"chr": "17", "start": 41196312, "end": 41277500},
    "TP53": {"chr": "17", "start": 7571720, "end": 7590868},
    "CFTR": {"chr": "7", "start": 117120017, "end": 117308719},
    "PKD1": {"chr": "16", "start": 2138711, "end": 2185899},
    "LDLR": {"chr": "19", "start": 11200038, "end": 11244505},
    "APOB": {"chr": "2", "start": 21224301, "end": 21266945},
    "BRCA2": {"chr": "13", "start": 32889617, "end": 32973809},
}

# VCF header
VCF_HEADER = """##fileformat=VCFv4.2
##fileDate={}
##source=create_test_vcf_for_regression.py
##reference=GRCh37
##contig=<ID=2,length=243199373>
##contig=<ID=7,length=159138663>
##contig=<ID=13,length=115169878>
##contig=<ID=16,length=90354753>
##contig=<ID=17,length=81195210>
##contig=<ID=19,length=59128983>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
##FORMAT=<ID=AD,Number=R,Type=Integer,Description="Allelic depths">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality">
##FORMAT=<ID=PL,Number=G,Type=Integer,Description="Phred-scaled genotype likelihoods">
##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">
##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles">
##INFO=<ID=DP,Number=1,Type=Integer,Description="Total depth">
##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations from SnpEff">
##INFO=<ID=LOF,Number=.,Type=String,Description="Loss of function prediction">
##INFO=<ID=CADD_phred,Number=A,Type=Float,Description="CADD Phred score">
##INFO=<ID=dbNSFP_gnomAD_exomes_AF,Number=A,Type=Float,Description="gnomAD exomes allele frequency">
##INFO=<ID=ClinVar_CLNSIG,Number=.,Type=String,Description="ClinVar clinical significance">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	SAMPLE001	SAMPLE002	SAMPLE003
""".format(datetime.now().strftime("%Y%m%d"))

# Test variants - mix of pathogenic, benign, rare, etc.
TEST_VARIANTS = [
    # BRCA1 variants
    {
        "chr": "17", "pos": 41245466, "ref": "G", "alt": "A",
        "qual": 2500, "info": "AC=1;AF=0.167;AN=6;DP=200;ANN=A|stop_gained|HIGH|BRCA1|BRCA1|transcript|NM_007294.3|protein_coding|11/24|c.1687C>T|p.Gln563*|1806/7207|1687/5655|563/1884||;LOF=(BRCA1|BRCA1|1|1.00);CADD_phred=42;dbNSFP_gnomAD_exomes_AF=0.00001;ClinVar_CLNSIG=Pathogenic",
        "samples": ["0/1:100,98:198:99:2500,0,2450", "0/0:120,0:120:99:0,360,3600", "0/0:110,0:110:99:0,330,3300"]
    },
    {
        "chr": "17", "pos": 41244936, "ref": "T", "alt": "C",
        "qual": 1800, "info": "AC=2;AF=0.333;AN=6;DP=180;ANN=C|missense_variant|MODERATE|BRCA1|BRCA1|transcript|NM_007294.3|protein_coding|11/24|c.2217A>G|p.Met739Val|2336/7207|2217/5655|739/1884||;CADD_phred=25;dbNSFP_gnomAD_exomes_AF=0.0005;ClinVar_CLNSIG=Uncertain_significance",
        "samples": ["0/1:90,88:178:99:1800,0,1750", "0/1:95,93:188:99:1850,0,1800", "0/0:100,0:100:99:0,300,3000"]
    },
    
    # TP53 variants
    {
        "chr": "17", "pos": 7577120, "ref": "C", "alt": "T",
        "qual": 3000, "info": "AC=1;AF=0.167;AN=6;DP=220;ANN=T|missense_variant|MODERATE|TP53|TP53|transcript|NM_000546.5|protein_coding|7/11|c.743G>A|p.Arg248Gln|885/2591|743/1182|248/393||;CADD_phred=35;dbNSFP_gnomAD_exomes_AF=0.000002;ClinVar_CLNSIG=Pathogenic",
        "samples": ["0/0:130,0:130:99:0,390,3900", "0/1:110,108:218:99:3000,0,2950", "0/0:120,0:120:99:0,360,3600"]
    },
    
    # CFTR variants
    {
        "chr": "7", "pos": 117199646, "ref": "A", "alt": "G",
        "qual": 2800, "info": "AC=2;AF=0.333;AN=6;DP=210;ANN=G|missense_variant|MODERATE|CFTR|CFTR|transcript|NM_000492.3|protein_coding|13/27|c.1652G>A|p.Gly551Asp|1784/6129|1652/4443|551/1480||;CADD_phred=28;dbNSFP_gnomAD_exomes_AF=0.0001;ClinVar_CLNSIG=Pathogenic",
        "samples": ["0/1:105,103:208:99:2800,0,2750", "0/0:115,0:115:99:0,345,3450", "0/1:100,98:198:99:2700,0,2650"]
    },
    
    # PKD1 variants
    {
        "chr": "16", "pos": 2160000, "ref": "G", "alt": "C",
        "qual": 2200, "info": "AC=1;AF=0.167;AN=6;DP=190;ANN=C|stop_gained|HIGH|PKD1|PKD1|transcript|NM_001009944.2|protein_coding|15/46|c.3850C>G|p.Arg1284*|4176/14149|3850/12915|1284/4304||;LOF=(PKD1|PKD1|1|1.00);CADD_phred=38;dbNSFP_gnomAD_exomes_AF=0.00005;ClinVar_CLNSIG=Likely_pathogenic",
        "samples": ["0/1:95,93:188:99:2200,0,2150", "0/0:105,0:105:99:0,315,3150", "0/0:100,0:100:99:0,300,3000"]
    },
    
    # LDLR variants
    {
        "chr": "19", "pos": 11216207, "ref": "C", "alt": "T",
        "qual": 2600, "info": "AC=1;AF=0.167;AN=6;DP=200;ANN=T|missense_variant|MODERATE|LDLR|LDLR|transcript|NM_000527.4|protein_coding|4/18|c.502G>A|p.Asp168Asn|681/5253|502/2586|168/861||;CADD_phred=24;dbNSFP_gnomAD_exomes_AF=0.001;ClinVar_CLNSIG=Pathogenic",
        "samples": ["0/0:110,0:110:99:0,330,3300", "0/1:100,98:198:99:2600,0,2550", "0/0:105,0:105:99:0,315,3150"]
    },
    
    # APOB variants
    {
        "chr": "2", "pos": 21232200, "ref": "G", "alt": "A",
        "qual": 2400, "info": "AC=2;AF=0.333;AN=6;DP=195;ANN=A|missense_variant|MODERATE|APOB|APOB|transcript|NM_000384.2|protein_coding|26/29|c.10580C>T|p.Arg3527Trp|10713/14121|10580/13716|3527/4571||;CADD_phred=26;dbNSFP_gnomAD_exomes_AF=0.0002;ClinVar_CLNSIG=Pathogenic",
        "samples": ["0/1:98,96:194:99:2400,0,2350", "0/1:97,95:192:99:2350,0,2300", "0/0:100,0:100:99:0,300,3000"]
    },
    
    # Additional BRCA2 variant (already in genes)
    {
        "chr": "13", "pos": 32914437, "ref": "AG", "alt": "A",
        "qual": 2900, "info": "AC=1;AF=0.167;AN=6;DP=215;ANN=A|frameshift_variant|HIGH|BRCA2|BRCA2|transcript|NM_000059.3|protein_coding|11/27|c.5946delT|p.Ser1982fs|6146/11388|5946/10257|1982/3418||;LOF=(BRCA2|BRCA2|1|1.00);CADD_phred=36;ClinVar_CLNSIG=Pathogenic",
        "samples": ["0/0:115,0:115:99:0,345,3450", "0/0:110,0:110:99:0,330,3300", "0/1:108,106:214:99:2900,0,2850"]
    },
    
    # Add some benign/common variants
    {
        "chr": "17", "pos": 41251800, "ref": "A", "alt": "G",
        "qual": 1500, "info": "AC=3;AF=0.5;AN=6;DP=150;ANN=G|synonymous_variant|LOW|BRCA1|BRCA1|transcript|NM_007294.3|protein_coding|9/24|c.1281T>C|p.Ser427Ser|1400/7207|1281/5655|427/1884||;CADD_phred=5;dbNSFP_gnomAD_exomes_AF=0.25;ClinVar_CLNSIG=Benign",
        "samples": ["0/1:75,73:148:99:1500,0,1450", "1/1:0,80:80:99:2400,240,0", "0/0:85,0:85:99:0,255,2550"]
    },
]


def main():
    """Create test VCF file with regression test variants."""
    import subprocess
    import os
    
    temp_vcf = "tests/fixtures/test_regression_variants.GRCh37.annotated.vcf"
    output_file = "tests/fixtures/test_regression_variants.GRCh37.annotated.vcf.gz"
    
    print(f"Creating test VCF file: {output_file}")
    
    # Sort variants by chromosome and position
    sorted_variants = sorted(TEST_VARIANTS, key=lambda v: (
        int(v["chr"]) if v["chr"].isdigit() else 99,
        int(v["pos"])
    ))
    
    # Write uncompressed VCF first
    with open(temp_vcf, "w") as f:
        # Write header
        f.write(VCF_HEADER)
        
        # Write sorted variants
        for variant in sorted_variants:
            line = "{}\t{}\t.\t{}\t{}\t{}\tPASS\t{}\tGT:AD:DP:GQ:PL\t{}\n".format(
                variant["chr"],
                variant["pos"],
                variant["ref"],
                variant["alt"],
                variant["qual"],
                variant["info"],
                "\t".join(variant["samples"])
            )
            f.write(line)
    
    # Compress with bgzip
    subprocess.run(["bgzip", "-f", temp_vcf], check=True)
    
    # Create tabix index
    subprocess.run(["tabix", "-p", "vcf", output_file], check=True)
    
    print(f"Created test VCF with {len(TEST_VARIANTS)} variants")
    print(f"Genes included: {', '.join(TEST_GENES.keys())}")
    
    # Also create the test gene file
    with open("tests/fixtures/test_genes.txt", "w") as f:
        f.write("BRCA1\n")
        f.write("BRCA2\n")
        f.write("TP53\n")
    print("Created test_genes.txt")
    
    # Create a simple PED file
    with open("tests/fixtures/test_family.ped", "w") as f:
        f.write("FAM001\tSAMPLE001\t0\t0\t1\t2\n")  # Affected child
        f.write("FAM001\tSAMPLE002\t0\t0\t2\t1\n")  # Unaffected mother
        f.write("FAM001\tSAMPLE003\t0\t0\t1\t1\n")  # Unaffected father
    print("Created test_family.ped")
    
    # Create a simple BED file
    with open("tests/fixtures/test_regions.bed", "w") as f:
        f.write("17\t41200000\t41300000\tBRCA1_region\n")
        f.write("13\t32800000\t33000000\tBRCA2_region\n")
    print("Created test_regions.bed")


if __name__ == "__main__":
    main()