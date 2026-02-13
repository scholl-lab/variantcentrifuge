"""Mock implementations of external bioinformatics tools."""


class MockBCFTools:
    """Mock bcftools for testing without dependencies."""

    def __init__(self):
        self.commands_executed = []

    def run(self, args: list[str], output_file: str | None = None) -> str:
        """Mock bcftools execution."""
        self.commands_executed.append(args)

        if args[0] == "view" and "-h" in args:
            # Mock header output
            return self._mock_vcf_header()
        elif args[0] == "view" and ("-R" in args or "-r" in args):
            # Mock variant extraction
            return self._mock_extracted_variants()
        elif args[0] == "query":
            # Mock query output
            return self._mock_query_output(args)
        else:
            return ""

    def _mock_vcf_header(self) -> str:
        """Return mock VCF header."""
        return """##fileformat=VCFv4.2
##contig=<ID=chr1,length=249250621>
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\tSample3
"""

    def _mock_extracted_variants(self) -> str:
        """Return mock extracted variants."""
        header = self._mock_vcf_header()
        variants = [
            "chr1\t100\t.\tA\tT\t100\tPASS\t.\tGT\t0/1\t0/0\t0/1",
            "chr1\t200\t.\tG\tC\t100\tPASS\t.\tGT\t0/0\t0/1\t1/1",
        ]
        return header + "\n".join(variants)

    def _mock_query_output(self, args: list[str]) -> str:
        """Return mock query output."""
        if "-l" in args:  # List samples
            return "Sample1\nSample2\nSample3\n"
        return ""


class MockSnpEff:
    """Mock snpEff for testing."""

    def __init__(self):
        self.commands_executed = []

    def run(self, args: list[str], output_file: str | None = None) -> str:
        """Mock snpEff execution."""
        self.commands_executed.append(args)

        if "genes2bed" in args:
            return self._mock_genes2bed(args)
        elif "ann" in args or "eff" in args:
            return self._mock_annotation()
        else:
            return ""

    def _mock_genes2bed(self, args: list[str]) -> str:
        """Return mock BED output for genes."""
        # Find the gene names in args
        genes = []
        for i, arg in enumerate(args):
            if arg == "-g" and i + 1 < len(args):
                genes = args[i + 1].split()
                break

        bed_lines = []
        for i, gene in enumerate(genes):
            # Mock BED coordinates
            start = 1000 * (i + 1)
            end = start + 500
            bed_lines.append(f"chr1\t{start}\t{end}\t{gene}\t0\t+")

        return "\n".join(bed_lines)

    def _mock_annotation(self) -> str:
        """Return mock annotated VCF."""
        return """##fileformat=VCFv4.2
##INFO=<ID=ANN,Number=.,Type=String,Description="Functional annotations">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tT\t100\tPASS\tANN=T|missense_variant|MODERATE|GENE1
"""


class MockSnpSift:
    """Mock SnpSift for testing."""

    def __init__(self):
        self.commands_executed = []

    def run(self, args: list[str], output_file: str | None = None) -> str:
        """Mock SnpSift execution."""
        self.commands_executed.append(args)

        if "filter" in args:
            return self._mock_filter()
        elif "extractFields" in args:
            return self._mock_extract_fields(args)
        elif "split" in args:
            return self._mock_split()
        else:
            return ""

    def _mock_filter(self) -> str:
        """Return mock filtered VCF."""
        return """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tT\t100\tPASS\tAF=0.01
"""

    def _mock_extract_fields(self, args: list[str]) -> str:
        """Return mock extracted fields TSV."""
        # Find fields in args
        fields = []
        for i, arg in enumerate(args):
            if arg == "-s" and i + 1 < len(args):
                fields = args[i + 1].split(",")
                break

        if not fields:
            # Default fields
            fields = ["CHROM", "POS", "REF", "ALT"]

        # Mock header
        header = "\t".join(fields)

        # Mock data rows
        rows = [
            "chr1\t100\tA\tT",
            "chr1\t200\tG\tC",
        ]

        return header + "\n" + "\n".join(rows)

    def _mock_split(self) -> str:
        """Return mock split VCF."""
        return """##fileformat=VCFv4.2
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tT\t100\tPASS\t.
chr1\t100\t.\tA\tG\t100\tPASS\t.
"""


class MockBedTools:
    """Mock bedtools for testing."""

    def __init__(self):
        self.commands_executed = []

    def run(self, args: list[str], output_file: str | None = None) -> str:
        """Mock bedtools execution."""
        self.commands_executed.append(args)

        if "intersect" in args:
            return self._mock_intersect()
        elif "merge" in args:
            return self._mock_merge()
        else:
            return ""

    def _mock_intersect(self) -> str:
        """Return mock intersection result."""
        return "chr1\t100\t200\tregion1\n"

    def _mock_merge(self) -> str:
        """Return mock merged BED."""
        return "chr1\t100\t500\n"
