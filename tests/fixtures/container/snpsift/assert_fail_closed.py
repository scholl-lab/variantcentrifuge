"""Exercise production SnpSift parser-diagnostic fail-closed behavior."""

from pathlib import Path
from tempfile import TemporaryDirectory

from variantcentrifuge.filters import apply_snpsift_filter

fixture = Path(__file__).with_name("input.vcf")
with TemporaryDirectory(prefix="vc-snpsift-") as temporary_directory:
    output = Path(temporary_directory) / "filtered.vcf.gz"
    try:
        apply_snpsift_filter(
            str(fixture),
            "(( QUAL >= 20 )",
            {"threads": 1},
            str(output),
        )
    except RuntimeError as error:
        expected_diagnostic = "SnpSift filter reported parser diagnostics"
        if expected_diagnostic not in str(error):
            raise RuntimeError("SnpSift failed without the required parser diagnostic") from error
    else:
        raise RuntimeError("SnpSift parser diagnostics did not fail closed")

print("SnpSift parser diagnostics failed closed")
