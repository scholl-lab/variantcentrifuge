"""Regression tests for #154: dead code removal and weights.py line budget."""

from pathlib import Path

import variantcentrifuge.display_utils
import variantcentrifuge.helpers
import variantcentrifuge.utils
from variantcentrifuge.association import weights


def test_helpers_public_exports_preserved():
    """Verify public exports in helpers.__all__ remain intact (#154)."""
    public_exports = [
        "annotate_variants_with_gene_lists",
        "assign_case_control_counts",
        "build_sample_phenotype_map",
        "check_file",
        "compute_phenotype_based_case_control_assignment",
        "determine_case_control_sets",
        "dump_df_to_xlsx",
        "extract_gencode_id",
        "extract_sample_and_genotype",
        "genotype_to_allele_count",
        "get_vcf_names",
        "get_vcf_regions",
        "get_vcf_samples",
        "get_vcf_size",
        "load_gene_list",
        "match_igv_link_columns",
        "read_sequencing_manifest",
    ]
    for func_name in public_exports:
        assert hasattr(variantcentrifuge.helpers, func_name), (
            f"{func_name} must remain importable in helpers.py"
        )
        assert func_name in variantcentrifuge.helpers.__all__, (
            f"{func_name} must remain in helpers.__all__"
        )


def test_utils_public_interface_preserved():
    """Verify documented public functions in utils remain importable (#154)."""
    assert hasattr(variantcentrifuge.utils, "check_external_tools"), (
        "check_external_tools must remain in utils.py as documented in api/utils.md"
    )
    assert variantcentrifuge.utils.check_external_tools([]) is True


def test_utils_retained_functions_present():
    """Verify live functions in utils are retained."""
    retained = [
        "compute_base_name",
        "generate_igv_safe_filename_base",
        "get_tool_version",
        "normalize_vcf_headers",
        "run_command",
        "smart_open",
    ]
    for func_name in retained:
        assert hasattr(variantcentrifuge.utils, func_name), f"{func_name} must remain in utils.py"


def test_display_utils_dead_functions_removed():
    """Verify 5 unused display functions are removed from display_utils."""
    removed_functions = [
        "display_dependency_graph",
        "display_stage_categories",
        "display_resume_safety_check",
        "format_stage_list_simple",
        "display_command_suggestions",
    ]
    for func_name in removed_functions:
        assert not hasattr(variantcentrifuge.display_utils, func_name), (
            f"{func_name} should be removed from display_utils.py"
        )


def test_display_utils_retained_functions_present():
    """Verify live display functions remain in display_utils."""
    retained = [
        "display_available_stages",
        "display_enhanced_status",
    ]
    for func_name in retained:
        assert hasattr(variantcentrifuge.display_utils, func_name), (
            f"{func_name} must remain in display_utils.py"
        )


def test_weights_physical_line_count():
    """Verify variantcentrifuge/association/weights.py is <= 600 physical lines."""
    weights_path = Path(weights.__file__)
    line_count = len(weights_path.read_text(encoding="utf-8").splitlines())
    assert line_count <= 600, (
        f"weights.py has {line_count} lines, which exceeds the 600-line budget"
    )


def test_weights_public_interface_preserved():
    """Verify all public functions in weights.py are preserved."""
    public_functions = [
        "beta_maf_weights",
        "combined_weights",
        "cadd_weights",
        "get_weights",
        "resolve_score_weight_column",
        "revel_weights",
        "score_column_weights",
        "uniform_weights",
    ]
    for func_name in public_functions:
        assert hasattr(weights, func_name), f"{func_name} must remain accessible in weights.py"


def test_loc_allowlist_retired_entries():
    """Verify Wave 1 and weights.py entries are retired from .loc-allowlist."""
    allowlist_path = Path(__file__).resolve().parent.parent.parent / ".loc-allowlist"
    assert allowlist_path.exists(), ".loc-allowlist must exist"

    lines = [
        line.strip()
        for line in allowlist_path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.strip().startswith("#")
    ]
    retired_files = [
        "variantcentrifuge/annotator.py",
        "variantcentrifuge/association/backends/python_backend.py",
        "variantcentrifuge/association/backends/r_backend.py",
        "variantcentrifuge/association/engine.py",
        "variantcentrifuge/association/weights.py",
        "variantcentrifuge/inheritance/vectorized_deducer.py",
        "variantcentrifuge/stages/processing_stages.py",
    ]
    allowlist_map = {line.split(":")[0]: int(line.split(":")[1]) for line in lines}

    for retired in retired_files:
        assert retired not in allowlist_map, f"{retired} should be retired from .loc-allowlist"

    # At most 6 Track 2A/2B files should remain (and eventually 0 when tracks merge)
    assert len(allowlist_map) <= 6, (
        f"Expected at most 6 entries in .loc-allowlist, got {len(allowlist_map)}: "
        f"{list(allowlist_map.keys())}"
    )
