"""Unit tests for argument quoting in Snakemake workflow."""

import re
import shlex
import subprocess
from pathlib import Path


def test_smk_rule_quotes_parameters() -> None:
    """Verify that workflow/rules/variantcentrifuge.smk uses :q and mkdir -p --."""
    smk_path = Path(__file__).resolve().parents[2] / "workflow" / "rules" / "variantcentrifuge.smk"
    content = smk_path.read_text()

    assert "mkdir -p -- {params.output_dir:q}" in content
    assert "-c {params.config_file:q}" in content
    assert "-g {params.genes:q}" in content
    assert "-v {input.vcf:q}" in content
    assert "--output-file {output.tsv:q}" in content
    assert "--output-dir {params.output_dir:q}" in content
    assert "--log-level {params.log_level:q}" in content
    assert "> {log:q} 2>&1" in content
    # extra_flags must NOT be quoted as a single argument
    assert "{params.extra_flags}" in content
    assert "{params.extra_flags:q}" not in content


def test_workflow_command_rendering_prevents_injection_and_preserves_boundaries(
    tmp_path: Path,
) -> None:
    """Simulate Snakemake :q rendering and execute with shell dummy wrapper."""
    smk_path = Path(__file__).resolve().parents[2] / "workflow" / "rules" / "variantcentrifuge.smk"
    content = smk_path.read_text()

    # Extract shell section
    shell_match = re.search(r'shell:\s*"""(.*?)"""', content, re.DOTALL)
    assert shell_match, "Shell block not found in variantcentrifuge.smk"
    shell_template = shell_match.group(1)

    marker_file = tmp_path / "injected.marker"
    params = {
        "params.output_dir": str(tmp_path / "out dir with spaces"),
        "params.config_file": "/path/to/config.json",
        "params.genes": f"PKD1$(touch {marker_file}) PKD2",
        "input.vcf": "/path with spaces/input.vcf.gz",
        "output.tsv": str(tmp_path / "out dir with spaces" / "sample.vc_analysis.tsv"),
        "params.log_level": "INFO",
        "threads": "4",
        "params.extra_flags": "--preset rare --append-extra-sample-fields 'GEN[*].DP'",
        "log": str(tmp_path / "out dir with spaces" / "log.txt"),
    }

    # Simulate Snakemake formatting where {field:q} -> shlex.quote(val) and {field} -> val
    rendered = shell_template
    for key, val in params.items():
        rendered = rendered.replace(f"{{{key}:q}}", shlex.quote(str(val)))
        rendered = rendered.replace(f"{{{key}}}", str(val))

    # Test executing in bash with dummy mkdir and variantcentrifuge functions
    # Dummy variantcentrifuge appends its argv to a log file
    args_log = tmp_path / "args.log"
    preamble = f"""
    mkdir() {{ /bin/mkdir -p "$@"; }}
    variantcentrifuge() {{
        for arg in "$@"; do
            echo "ARG:$arg" >> {shlex.quote(str(args_log))}
        done
    }}
    """
    script = preamble + "\n" + rendered

    res = subprocess.run(["bash", "-c", script], capture_output=True, text=True)
    assert res.returncode == 0, f"Bash script failed: {res.stderr}"
    assert not marker_file.exists(), "Command substitution was executed in workflow shell!"

    args = args_log.read_text().splitlines()
    # Check that -g received the entire multi-gene string as one argument
    g_idx = args.index("ARG:-g")
    assert args[g_idx + 1] == f"ARG:PKD1$(touch {marker_file}) PKD2"

    # Check extra_flags preserved boundaries
    assert "ARG:--preset" in args
    assert "ARG:rare" in args
    assert "ARG:--append-extra-sample-fields" in args
    assert "ARG:GEN[*].DP" in args
