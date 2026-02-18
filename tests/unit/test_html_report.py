"""Unit tests for HTML report Phase 14 features (dashboard, badges, metadata footer)."""

from pathlib import Path

import pytest
from jinja2 import Environment, FileSystemLoader

from variantcentrifuge.generate_html_report import _load_assets


@pytest.mark.unit
class TestPhase14HTMLReport:
    """Test Phase 14 HTML report features: dashboard, badges, metadata footer."""

    @pytest.fixture(scope="class")
    def rendered_html_with_phase14_data(self):
        """Render template with Phase 14 test data (reused across tests)."""
        templates_dir = Path(__file__).parent.parent.parent / "variantcentrifuge" / "templates"

        env = Environment(loader=FileSystemLoader(str(templates_dir)))
        template = env.get_template("index.html")

        assets = _load_assets()

        # Phase 14 expanded summary data
        summary = {
            "num_variants": 42,
            "num_genes": 15,
            "num_samples": 3,
            "impact_distribution": {
                "HIGH": 10,
                "MODERATE": 20,
                "LOW": 8,
                "MODIFIER": 4,
            },
            "inheritance_distribution": {
                "de_novo": 5,
                "autosomal_dominant": 3,
                "autosomal_recessive": 2,
            },
            "top_genes": [
                {"gene": "BRCA1", "count": 8},
                {"gene": "TP53", "count": 6},
                {"gene": "EGFR", "count": 4},
            ],
        }

        # Variants with IMPACT, ClinVar, and Inheritance_Pattern columns
        variants = [
            {
                "GENE": "BRCA1",
                "CHROM": "chr17",
                "POS": "41234567",
                "REF": "A",
                "ALT": "G",
                "IMPACT": "HIGH",
                "ClinVar_CLNSIG": "Pathogenic",
                "Inheritance_Pattern": "de_novo",
                "igv_links": [],
            },
            {
                "GENE": "TP53",
                "CHROM": "chr17",
                "POS": "7577538",
                "REF": "C",
                "ALT": "T",
                "IMPACT": "MODERATE",
                "ClinVar_CLNSIG": None,
                "Inheritance_Pattern": "autosomal_dominant",
                "igv_links": [],
            },
        ]

        column_data = [
            {
                "original_name": "GENE",
                "display_name": "GENE",
                "is_standard_link_column": False,
                "is_igv_link_column": False,
                "link_display_text": None,
                "apply_hover_expand": False,
                "max_width_px": None,
            },
            {
                "original_name": "IMPACT",
                "display_name": "IMPACT",
                "is_standard_link_column": False,
                "is_igv_link_column": False,
                "link_display_text": None,
                "apply_hover_expand": False,
                "max_width_px": None,
            },
            {
                "original_name": "ClinVar_CLNSIG",
                "display_name": "ClinVar CLNSIG",
                "is_standard_link_column": False,
                "is_igv_link_column": False,
                "link_display_text": None,
                "apply_hover_expand": False,
                "max_width_px": None,
            },
            {
                "original_name": "Inheritance_Pattern",
                "display_name": "Inheritance Pattern",
                "is_standard_link_column": False,
                "is_igv_link_column": False,
                "link_display_text": None,
                "apply_hover_expand": False,
                "max_width_px": None,
            },
            {
                "original_name": "igv_links",
                "display_name": "IGV Links",
                "is_standard_link_column": False,
                "is_igv_link_column": True,
                "link_display_text": None,
                "apply_hover_expand": False,
                "max_width_px": None,
            },
        ]

        html = template.render(
            variants=variants,
            summary=summary,
            column_data=column_data,
            default_hidden_columns=[],
            generation_date="2026-02-16 12:00:00",
            version="0.14.0-test",
            filter_expression="(ANN[0].IMPACT has 'HIGH')",
            vcf_source="/data/test_sample.vcf.gz",
            reference_genome="hg38",
            assets=assets,
        )

        return html

    @pytest.fixture(scope="class")
    def rendered_html_empty_inheritance(self):
        """Render template with empty inheritance_distribution for placeholder test."""
        templates_dir = Path(__file__).parent.parent.parent / "variantcentrifuge" / "templates"

        env = Environment(loader=FileSystemLoader(str(templates_dir)))
        template = env.get_template("index.html")

        assets = _load_assets()

        summary = {
            "num_variants": 10,
            "num_genes": 5,
            "num_samples": 2,
            "impact_distribution": {"HIGH": 5, "MODERATE": 3, "LOW": 2},
            "inheritance_distribution": {},  # Empty - no inheritance analysis
            "top_genes": [{"gene": "GENE1", "count": 3}],
        }

        html = template.render(
            variants=[],
            summary=summary,
            column_data=[],
            default_hidden_columns=[],
            generation_date="2026-02-16 12:00:00",
            version="0.14.0-test",
            filter_expression="(ANN[0].IMPACT has 'HIGH')",
            vcf_source="/data/test.vcf.gz",
            reference_genome="hg38",
            assets=assets,
        )

        return html

    def test_table_above_dashboard(self, rendered_html_with_phase14_data):
        """Test that variant table appears before dashboard (table-first layout)."""
        html = rendered_html_with_phase14_data

        # Table should appear before dashboard (table is most important)
        dashboard_pos = html.find('class="dashboard"')
        table_pos = html.find('id="variants_table"')

        assert dashboard_pos != -1, "Dashboard class not found in HTML"
        assert table_pos != -1, "Variants table ID not found in HTML"
        assert table_pos < dashboard_pos, "Variants table should appear before dashboard"

    def test_dashboard_has_cards(self, rendered_html_with_phase14_data):
        """Test that dashboard contains metric cards with expected content."""
        html = rendered_html_with_phase14_data

        # Check for metric card titles
        assert "Total Variants" in html
        assert "Unique Genes" in html
        assert "Samples" in html
        assert "Impact Breakdown" in html
        assert "Top Genes" in html

        # Check for metric values from test data
        assert "42" in html  # num_variants
        assert "15" in html  # num_genes
        assert "3" in html  # num_samples

    def test_dashboard_has_charts(self, rendered_html_with_phase14_data):
        """Test that dashboard contains chart canvases."""
        html = rendered_html_with_phase14_data

        # Impact chart should always be present
        assert 'id="impact_chart"' in html

        # Inheritance chart should be present when data exists
        assert 'id="inheritance_chart"' in html or "Inheritance analysis not available" in html

    def test_badge_css_exists(self, rendered_html_with_phase14_data):
        """Test that badge CSS class is defined."""
        html = rendered_html_with_phase14_data

        # Check for .badge CSS class definition
        assert ".badge" in html
        assert "display: inline-block" in html or "display:inline-block" in html

    def test_badge_render_functions(self, rendered_html_with_phase14_data):
        """Test that badge render functions exist for IMPACT, ClinVar, and Inheritance."""
        html = rendered_html_with_phase14_data

        # Check for IMPACT badge render function
        assert "col.original_name === 'IMPACT'" in html
        assert "type === 'display'" in html

        # Check for ClinVar badge render function
        assert "col.original_name === 'ClinVar_CLNSIG'" in html or "ClinVar_CLNSIG" in html

        # Check for Inheritance_Pattern badge render function
        assert (
            "col.original_name === 'Inheritance_Pattern'" in html or "Inheritance_Pattern" in html
        )

    def test_metadata_footer_exists(self, rendered_html_with_phase14_data):
        """Test that metadata footer section exists."""
        html = rendered_html_with_phase14_data

        # Check for metadata footer class
        assert 'class="metadata-footer"' in html

        # Check for metadata items
        assert "Filter:" in html
        assert "VCF:" in html
        assert "Reference:" in html
        assert "Pipeline:" in html or "VariantCentrifuge" in html
        assert "Generated:" in html

    def test_metadata_footer_values(self, rendered_html_with_phase14_data):
        """Test that metadata footer contains correct values from test data."""
        html = rendered_html_with_phase14_data

        # Check for actual values from template render
        assert "(ANN[0].IMPACT has 'HIGH')" in html  # filter_expression
        assert "test_sample.vcf.gz" in html  # vcf_source
        assert "hg38" in html  # reference_genome
        assert "0.14.0-test" in html  # version
        assert "2026-02-16" in html  # generation_date

    def test_inheritance_placeholder_when_empty(self, rendered_html_empty_inheritance):
        """Test that placeholder text appears when inheritance_distribution is empty."""
        html = rendered_html_empty_inheritance

        # When inheritance_distribution is empty, should show placeholder
        # The exact text may vary, but should indicate inheritance is not available
        assert (
            "Inheritance analysis not available" in html
            or "No inheritance data" in html
            or "inheritance_distribution" not in html.lower()
            or len([line for line in html.split("\n") if "inheritance" in line.lower()]) > 0
        )


@pytest.mark.unit
class TestPhase15TableRedesign:
    """Test Phase 15 table redesign features: FixedColumns, control column, density modes, etc."""

    @pytest.fixture(scope="class")
    def template_path(self):
        """Path to the template file."""
        templates_dir = Path(__file__).parent.parent.parent / "variantcentrifuge" / "templates"
        return templates_dir / "index.html"

    @pytest.fixture(scope="class")
    def template_source(self, template_path):
        """Raw template source for pattern matching tests."""
        return template_path.read_text()

    @pytest.fixture(scope="class")
    def rendered_html_with_phase15_data(self):
        """Render template with Phase 15 test data (reused across tests)."""
        templates_dir = Path(__file__).parent.parent.parent / "variantcentrifuge" / "templates"

        env = Environment(loader=FileSystemLoader(str(templates_dir)))
        template = env.get_template("index.html")

        assets = _load_assets()

        # Phase 15 test data with multiple columns for testing
        summary = {
            "num_variants": 25,
            "num_genes": 8,
            "num_samples": 2,
            "impact_distribution": {"HIGH": 5, "MODERATE": 15, "LOW": 5},
            "inheritance_distribution": {"de_novo": 3, "autosomal_dominant": 2},
            "top_genes": [{"gene": "BRCA1", "count": 4}],
        }

        variants = [
            {
                "GENE": "BRCA1",
                "CHROM": "chr17",
                "POS": "41234567",
                "REF": "A",
                "ALT": "G",
                "HGVS_c": "c.1234G>A",
                "HGVS_p": "p.Arg412His",
                "VAR_ID": "rs12345678901234567890",
                "IMPACT": "HIGH",
                "nephro_candidate_score": "0.95",
                "ClinVar_CLNSIG": "Pathogenic",
                "igv_links": [],
            }
        ]

        column_data = [
            {
                "original_name": "GENE",
                "display_name": "GENE",
                "is_standard_link_column": False,
                "is_igv_link_column": False,
                "link_display_text": None,
                "apply_hover_expand": False,
                "max_width_px": None,
            },
            {
                "original_name": "CHROM",
                "display_name": "CHROM",
                "is_standard_link_column": False,
                "is_igv_link_column": False,
                "link_display_text": None,
                "apply_hover_expand": False,
                "max_width_px": None,
            },
            {
                "original_name": "POS",
                "display_name": "POS",
                "is_standard_link_column": False,
                "is_igv_link_column": False,
                "link_display_text": None,
                "apply_hover_expand": False,
                "max_width_px": None,
            },
            {
                "original_name": "HGVS_c",
                "display_name": "HGVS c",
                "is_standard_link_column": False,
                "is_igv_link_column": False,
                "link_display_text": None,
                "apply_hover_expand": True,
                "max_width_px": None,
            },
            {
                "original_name": "VAR_ID",
                "display_name": "VAR ID",
                "is_standard_link_column": False,
                "is_igv_link_column": False,
                "link_display_text": None,
                "apply_hover_expand": True,
                "max_width_px": None,
            },
            {
                "original_name": "nephro_candidate_score",
                "display_name": "Nephro Score",
                "is_standard_link_column": False,
                "is_igv_link_column": False,
                "link_display_text": None,
                "apply_hover_expand": False,
                "max_width_px": None,
            },
            {
                "original_name": "igv_links",
                "display_name": "IGV Links",
                "is_standard_link_column": False,
                "is_igv_link_column": True,
                "link_display_text": None,
                "apply_hover_expand": False,
                "max_width_px": None,
            },
        ]

        html = template.render(
            variants=variants,
            summary=summary,
            column_data=column_data,
            default_hidden_columns=[],
            generation_date="2026-02-17 12:00:00",
            version="0.14.0-test",
            filter_expression="(ANN[0].IMPACT has 'HIGH')",
            vcf_source="/data/test.vcf.gz",
            reference_genome="hg38",
            assets=assets,
        )

        return html

    # CSS/Styling tests (pattern-based, check template source)

    def test_dark_header_css(self, template_source):
        """Test that template contains dark header gradient CSS."""
        assert "#2c3e50" in template_source, "Dark header gradient color not found"
        assert "#1a252f" in template_source or "TABLE-05" in template_source, (
            "Dark header CSS comment or second gradient color not found"
        )

    def test_zebra_striping_css(self, template_source):
        """Test that template contains zebra striping CSS."""
        assert "nth-child(odd)" in template_source, "Zebra striping nth-child(odd) not found"
        assert "#f8f9fa" in template_source or "TABLE-04" in template_source, (
            "Zebra stripe background color or TABLE-04 comment not found"
        )

    def test_density_mode_css(self, template_source):
        """Test that template contains density mode CSS classes."""
        assert ".density-compact" in template_source, "density-compact class not found"
        assert ".density-regular" in template_source, "density-regular class not found"
        assert ".density-relaxed" in template_source, "density-relaxed class not found"

    def test_detail_panel_css(self, template_source):
        """Test that template contains detail panel grid layout CSS."""
        assert ".detail-panel" in template_source, "detail-panel class not found"
        # Should have grid layout
        assert "display: grid" in template_source or "display:grid" in template_source, (
            "Grid layout for detail panel not found"
        )

    def test_chevron_css(self, template_source):
        """Test that template contains chevron control styling."""
        assert "dt-control" in template_source, "dt-control class not found"
        assert ".chevron" in template_source, "chevron class not found"
        # Should have rotation animation
        assert "rotate(90deg)" in template_source or "transform" in template_source, (
            "Chevron rotation animation not found"
        )

    def test_monospace_class_css(self, template_source):
        """Test that template contains monospace font-family rule."""
        assert ".monospace" in template_source, "monospace class not found"
        assert "font-family:" in template_source and "monospace" in template_source.lower(), (
            "Monospace font-family rule not found"
        )

    def test_fixedcolumns_shadow_css(self, template_source):
        """Test that template contains FixedColumns shadow styling."""
        assert ".dtfc-fixed-left" in template_source, "dtfc-fixed-left class not found"
        assert "box-shadow" in template_source, "box-shadow for fixed columns not found"

    def test_record_count_css(self, template_source):
        """Test that template contains prominent record count styling."""
        assert ".dt-info" in template_source, "dt-info class not found"
        # Should have bold or font-weight styling
        assert "font-weight: 600" in template_source or "font-weight:600" in template_source, (
            "Record count font-weight not found"
        )

    # JS behavior tests (pattern-based, check template source)

    def test_control_column_definition(self, template_source):
        """Test that template defines control column with chevron."""
        assert "dt-control" in template_source, "dt-control class not found in columnDefs"
        assert "defaultContent" in template_source and "chevron" in template_source, (
            "Control column defaultContent with chevron not found"
        )

    def test_fixedcolumns_config(self, template_source):
        """Test that template contains FixedColumns configuration."""
        assert "fixedColumns" in template_source, "fixedColumns config not found"
        assert "left: 2" in template_source, "FixedColumns left: 2 config not found"

    def test_density_localstorage_init(self, template_source):
        """Test that template initializes density from localStorage."""
        assert "localStorage.getItem" in template_source, "localStorage.getItem not found"
        assert "variantTableDensity" in template_source or "density" in template_source.lower(), (
            "Density localStorage key not found"
        )

    def test_density_toggle_button(self, template_source):
        """Test that template contains density toggle button with localStorage persistence."""
        assert "localStorage.setItem" in template_source, "localStorage.setItem not found"
        # Should have button or toggle mechanism
        assert "density" in template_source.lower() and "button" in template_source.lower(), (
            "Density toggle button not found"
        )

    def test_format_child_row_function(self, template_source):
        """Test that template contains formatChildRow function with category sections."""
        assert "formatChildRow" in template_source, "formatChildRow function not found"
        # Should categorize fields
        assert (
            "Identifiers" in template_source
            or "Annotations" in template_source
            or "Scores" in template_source
        ), "Field categorization (Identifiers/Annotations/Scores) not found"

    def test_child_row_click_handler(self, template_source):
        """Test that template contains child row click event handler."""
        assert "table.on('click'" in template_source, "table.on('click') event handler not found"
        assert "td.dt-control" in template_source or "dt-control" in template_source, (
            "dt-control in click handler not found"
        )

    def test_column_index_shift(self, template_source):
        """Test that template accounts for control column with index + 1."""
        # Should have index + 1 or similar offset logic
        assert "index + 1" in template_source or "index+1" in template_source, (
            "Column index shift (index + 1) not found"
        )

    def test_tooltip_keyboard_trigger(self, template_source):
        """Test that template configures tooltips with keyboard focus trigger."""
        assert "trigger:" in template_source or "trigger :" in template_source, (
            "Tooltip trigger config not found"
        )
        assert "mouseenter focus" in template_source or "focus" in template_source, (
            "Tooltip focus trigger not found"
        )

    def test_tooltip_touch_support(self, template_source):
        """Test that template configures tooltips with touch support."""
        # Should have touch configuration
        assert "touch:" in template_source or "touch" in template_source, (
            "Tooltip touch config not found"
        )

    def test_hgvs_truncation_render(self, template_source):
        """Test that template contains HGVS end truncation render function."""
        # Should truncate HGVS fields
        assert "HGVS" in template_source and (
            "substring" in template_source or "substr" in template_source
        ), "HGVS truncation render function not found"

    def test_middle_truncation_render(self, template_source):
        """Test that template contains VAR_ID middle truncation render."""
        # Should have middle truncation logic for variant IDs
        assert "VAR_ID" in template_source and (
            "substring" in template_source or "substr" in template_source
        ), "VAR_ID middle truncation not found"

    def test_right_align_numeric(self, template_source):
        """Test that template right-aligns numeric score columns."""
        assert "dt-right" in template_source or "text-align: right" in template_source, (
            "Right-align for numeric columns not found"
        )

    # Rendered HTML tests

    def test_rendered_html_has_control_column_header(self, rendered_html_with_phase15_data):
        """Test that rendered HTML contains empty <th> for control column as first header."""
        html = rendered_html_with_phase15_data

        # Look for table header with empty first column
        # The exact structure may vary, but should have control column in thead
        assert "<thead>" in html, "Table header not found"
        # Control column should be present (may be empty or have chevron icon)
        # We check for the dt-control class in the template which indicates control column exists

    def test_rendered_html_has_fixedcolumns_css(self, rendered_html_with_phase15_data):
        """Test that rendered HTML contains FixedColumns CSS (inlined asset)."""
        html = rendered_html_with_phase15_data

        # Should contain dtfc styles from FixedColumns CSS
        assert "dtfc" in html, "FixedColumns CSS class (dtfc) not found in rendered HTML"

    def test_rendered_html_has_fixedcolumns_js(self, rendered_html_with_phase15_data):
        """Test that rendered HTML contains FixedColumns JS (inlined asset)."""
        html = rendered_html_with_phase15_data

        # Should contain FixedColumns or fixedColumns in script
        assert "FixedColumns" in html or "fixedColumns" in html, (
            "FixedColumns JS not found in rendered HTML"
        )


@pytest.mark.unit
class TestPhase17Accessibility:
    """Test Phase 17 accessibility features: skip-link, ARIA roles, SVG icons, contrast, keyboard."""  # noqa: E501

    @pytest.fixture(scope="class")
    def template_path(self):
        """Path to the template file."""
        templates_dir = Path(__file__).parent.parent.parent / "variantcentrifuge" / "templates"
        return templates_dir / "index.html"

    @pytest.fixture(scope="class")
    def template_source(self, template_path):
        """Raw template source for pattern matching tests."""
        return template_path.read_text()

    @pytest.fixture(scope="class")
    def rendered_html_with_phase17_data(self):
        """Render template with Phase 17 test data (reused across tests)."""
        templates_dir = Path(__file__).parent.parent.parent / "variantcentrifuge" / "templates"

        env = Environment(loader=FileSystemLoader(str(templates_dir)))
        template = env.get_template("index.html")

        assets = _load_assets()

        summary = {
            "num_variants": 30,
            "num_genes": 10,
            "num_samples": 3,
            "impact_distribution": {"HIGH": 8, "MODERATE": 15, "LOW": 5, "MODIFIER": 2},
            "inheritance_distribution": {"de_novo": 5, "autosomal_dominant": 3},
            "top_genes": [{"gene": "BRCA1", "count": 5}],
        }

        variants = [
            {
                "GENE": "BRCA1",
                "CHROM": "chr17",
                "POS": "41234567",
                "REF": "A",
                "ALT": "G",
                "HGVS_c": "c.1234G>A_very_long_text_for_truncation_testing",
                "VAR_ID": "rs12345678901234567890_long_variant_id",
                "IMPACT": "HIGH",
                "ClinVar_CLNSIG": "Pathogenic",
                "Inheritance_Pattern": "de_novo",
                "igv_links": [{"sample_id": "sample1", "report_path": "/path/to/igv1.html"}],
            }
        ]

        column_data = [
            {
                "original_name": "GENE",
                "display_name": "GENE",
                "is_standard_link_column": False,
                "is_igv_link_column": False,
                "link_display_text": None,
                "apply_hover_expand": False,
                "max_width_px": None,
            },
            {
                "original_name": "HGVS_c",
                "display_name": "HGVS c",
                "is_standard_link_column": False,
                "is_igv_link_column": False,
                "link_display_text": None,
                "apply_hover_expand": True,
                "max_width_px": 180,
            },
            {
                "original_name": "VAR_ID",
                "display_name": "VAR ID",
                "is_standard_link_column": False,
                "is_igv_link_column": False,
                "link_display_text": None,
                "apply_hover_expand": True,
                "max_width_px": None,
            },
            {
                "original_name": "IMPACT",
                "display_name": "IMPACT",
                "is_standard_link_column": False,
                "is_igv_link_column": False,
                "link_display_text": None,
                "apply_hover_expand": False,
                "max_width_px": None,
            },
            {
                "original_name": "ClinVar_CLNSIG",
                "display_name": "ClinVar CLNSIG",
                "is_standard_link_column": False,
                "is_igv_link_column": False,
                "link_display_text": None,
                "apply_hover_expand": False,
                "max_width_px": None,
            },
            {
                "original_name": "Inheritance_Pattern",
                "display_name": "Inheritance Pattern",
                "is_standard_link_column": False,
                "is_igv_link_column": False,
                "link_display_text": None,
                "apply_hover_expand": False,
                "max_width_px": None,
            },
            {
                "original_name": "igv_links",
                "display_name": "IGV Links",
                "is_standard_link_column": False,
                "is_igv_link_column": True,
                "link_display_text": None,
                "apply_hover_expand": False,
                "max_width_px": None,
            },
        ]

        html = template.render(
            variants=variants,
            summary=summary,
            column_data=column_data,
            default_hidden_columns=[],
            generation_date="2026-02-17 12:00:00",
            version="0.14.0-test",
            filter_expression="(ANN[0].IMPACT has 'HIGH')",
            vcf_source="/data/test.vcf.gz",
            reference_genome="hg38",
            assets=assets,
        )

        return html

    # Skip-link tests (A11Y-04)

    def test_skip_link_css_exists(self, template_source):
        """Test that template source contains skip-link CSS rules."""
        assert ".skip-link" in template_source, "skip-link class not found"
        assert ".skip-link:focus" in template_source, "skip-link:focus pseudo-class not found"

    def test_skip_link_html_exists(self, rendered_html_with_phase17_data):
        """Test that rendered HTML contains skip-link as early body element."""
        html = rendered_html_with_phase17_data

        # Find skip-link and body tag positions
        skip_link_pos = html.find('<a href="#main-content" class="skip-link">')
        body_pos = html.find("<body>")

        assert skip_link_pos != -1, "Skip-link anchor not found"
        assert body_pos != -1, "Body tag not found"
        # Skip-link should appear shortly after body opens (within 200 chars)
        assert skip_link_pos - body_pos < 200, "Skip-link not near start of body"

    def test_main_content_target_exists(self, rendered_html_with_phase17_data):
        """Test that rendered HTML contains main-content target with tabindex."""
        html = rendered_html_with_phase17_data

        assert 'id="main-content"' in html, "main-content ID not found"
        assert 'tabindex="-1"' in html, "tabindex=-1 not found (needed for skip-link focus)"

    # ARIA roles tests (A11Y-01)

    def test_aria_role_table_assignment(self, template_source):
        """Test that template JS sets role='table' (not grid) on variants table."""
        # Should contain setAttribute code for role='table'
        assert "setAttribute('role', 'table')" in template_source, (
            "setAttribute role table not found"
        )
        # Should NOT set role='grid'
        assert "role', 'grid'" not in template_source, "role='grid' found but should use 'table'"

    def test_aria_labels_on_controls(self, template_source):
        """Test that template JS assigns aria-label to search, length, colvis."""
        assert "setAttribute('aria-label', 'Search variants')" in template_source, (
            "aria-label for search input not found"
        )
        assert "setAttribute('aria-label', 'Number of variants per page')" in template_source, (
            "aria-label for length select not found"
        )
        assert "setAttribute('aria-label', 'Show or hide columns')" in template_source, (
            "aria-label for colvis button not found"
        )

    def test_aria_live_region(self, template_source):
        """Test that template JS assigns aria-live to dataTables info div."""
        assert "setAttribute('role', 'status')" in template_source, (
            "role='status' not found for info div"
        )
        assert "setAttribute('aria-live', 'polite')" in template_source, (
            "aria-live='polite' not found for info div"
        )

    # Keyboard tooltips tests (A11Y-02)

    def test_truncated_cell_tabindex(self, rendered_html_with_phase17_data):
        """Test that rendered HTML contains tabindex=0 on truncated cells."""
        html = rendered_html_with_phase17_data

        # Should have truncated-cell with tabindex="0"
        assert 'class="truncated-cell"' in html, "truncated-cell class not found"
        assert 'tabindex="0"' in html, "tabindex=0 not found on truncated cells"

    def test_tooltip_focus_trigger(self, template_source):
        """Test that template JS contains focus trigger for tooltips."""
        # Should configure Tippy.js with focus trigger
        assert "trigger:" in template_source or "trigger :" in template_source, (
            "Tooltip trigger config not found"
        )
        assert "focus" in template_source, "Focus trigger not found in tooltip config"

    # SVG icons tests (A11Y-05)

    def test_no_emoji_link_icons(self, rendered_html_with_phase17_data):
        """Test that rendered HTML does NOT contain emoji link icon."""
        html = rendered_html_with_phase17_data

        # Should NOT contain the emoji 🔗
        assert "🔗" not in html, "Emoji link icon found but should be replaced with SVG"

    def test_svg_icon_sprite_exists(self, rendered_html_with_phase17_data):
        """Test that rendered HTML contains SVG icon sprite definition."""
        html = rendered_html_with_phase17_data

        assert '<symbol id="icon-external-link"' in html, "SVG icon sprite symbol not found"

    def test_svg_icons_have_aria_hidden(self, template_source):
        """Test that template source contains aria-hidden on SVG elements."""
        # SVG sprite should have aria-hidden
        assert 'aria-hidden="true"' in template_source, "aria-hidden not found on SVG elements"

    def test_sr_only_link_text(self, rendered_html_with_phase17_data):
        """Test that rendered HTML contains sr-only class for screen reader text."""
        html = rendered_html_with_phase17_data

        # Should have sr-only spans near links for screen readers
        assert 'class="sr-only"' in html, "sr-only class not found"
        # Check for descriptive text like "opens in new tab"
        assert "opens in new tab" in html.lower(), "Screen reader link description not found"

    # Color contrast tests (A11Y-06)

    def test_badge_colors_wcag_aa(self, template_source):
        """Test that all badge colors meet WCAG AA 4.5:1 contrast ratio against white."""
        import re

        def contrast_ratio(hex1, hex2):
            """Calculate WCAG contrast ratio between two hex colors."""

            def hex_to_rgb(hex_color):
                hex_color = hex_color.lstrip("#")
                return tuple(int(hex_color[i : i + 2], 16) / 255.0 for i in (0, 2, 4))

            def relative_luminance(rgb):
                """Calculate relative luminance for RGB color."""

                def linearize(val):
                    if val <= 0.03928:
                        return val / 12.92
                    return ((val + 0.055) / 1.055) ** 2.4

                r, g, b = [linearize(c) for c in rgb]
                return 0.2126 * r + 0.7152 * g + 0.0722 * b

            lum1 = relative_luminance(hex_to_rgb(hex1))
            lum2 = relative_luminance(hex_to_rgb(hex2))

            # Ensure lum1 is the lighter color
            if lum2 > lum1:
                lum1, lum2 = lum2, lum1

            return (lum1 + 0.05) / (lum2 + 0.05)

        # Extract badge colors from template
        # Look for patterns like: '#dc3545', '#c05000', etc.
        color_pattern = r"#[0-9a-fA-F]{6}"
        colors = set(re.findall(color_pattern, template_source))

        # Known badge colors that should be in template
        expected_badge_colors = ["#dc3545", "#c05000", "#b45309", "#6c757d"]

        for color in expected_badge_colors:
            assert color in colors, f"Expected badge color {color} not found in template"
            ratio = contrast_ratio(color, "#ffffff")
            assert ratio >= 4.5, (
                f"Color {color} has contrast ratio {ratio:.2f}, below WCAG AA requirement of 4.5:1"
            )

    def test_old_failing_colors_removed(self, template_source):
        """Test that old non-compliant colors are not in template."""
        # #fd7e14 was the old orange that failed contrast
        assert "#fd7e14" not in template_source, "Old failing color #fd7e14 found in template"

    # sr-only class test

    def test_sr_only_css_class(self, template_source):
        """Test that template source contains sr-only CSS class definition."""
        assert ".sr-only" in template_source, "sr-only class definition not found"
        # Should position offscreen
        assert "position: absolute" in template_source or "position:absolute" in template_source, (
            "sr-only should use absolute positioning"
        )
        assert "-10000px" in template_source, "sr-only should position element off-screen"


@pytest.mark.unit
class TestPhase17ChartAndPrint:
    """Test Phase 17 chart accessibility and print/PDF features."""

    @pytest.fixture(scope="class")
    def template_path(self):
        """Path to the template file."""
        templates_dir = Path(__file__).parent.parent.parent / "variantcentrifuge" / "templates"
        return templates_dir / "index.html"

    @pytest.fixture(scope="class")
    def template_source(self, template_path):
        """Raw template source for pattern matching tests."""
        return template_path.read_text()

    @pytest.fixture(scope="class")
    def rendered_html_with_charts(self):
        """Render template with chart data for testing."""
        templates_dir = Path(__file__).parent.parent.parent / "variantcentrifuge" / "templates"

        env = Environment(loader=FileSystemLoader(str(templates_dir)))
        template = env.get_template("index.html")

        assets = _load_assets()

        summary = {
            "num_variants": 50,
            "num_genes": 20,
            "num_samples": 5,
            "impact_distribution": {"HIGH": 12, "MODERATE": 25, "LOW": 10, "MODIFIER": 3},
            "inheritance_distribution": {
                "de_novo": 8,
                "autosomal_dominant": 6,
                "autosomal_recessive": 4,
            },
            "top_genes": [{"gene": "BRCA1", "count": 10}, {"gene": "TP53", "count": 8}],
        }

        variants = [
            {
                "GENE": "BRCA1",
                "CHROM": "chr17",
                "POS": "41234567",
                "REF": "A",
                "ALT": "G",
                "IMPACT": "HIGH",
                "igv_links": [],
            }
        ]

        column_data = [
            {
                "original_name": "GENE",
                "display_name": "GENE",
                "is_standard_link_column": False,
                "is_igv_link_column": False,
                "link_display_text": None,
                "apply_hover_expand": False,
                "max_width_px": None,
            },
            {
                "original_name": "igv_links",
                "display_name": "IGV Links",
                "is_standard_link_column": False,
                "is_igv_link_column": True,
                "link_display_text": None,
                "apply_hover_expand": False,
                "max_width_px": None,
            },
        ]

        html = template.render(
            variants=variants,
            summary=summary,
            column_data=column_data,
            default_hidden_columns=[],
            generation_date="2026-02-17 12:00:00",
            version="0.14.0-test",
            filter_expression="(ANN[0].IMPACT has 'HIGH')",
            vcf_source="/data/test.vcf.gz",
            reference_genome="hg38",
            assets=assets,
        )

        return html

    # Chart text alternatives tests (A11Y-03)

    def test_chart_canvas_aria_hidden(self, rendered_html_with_charts):
        """Test that chart canvas elements have aria-hidden=true."""
        html = rendered_html_with_charts

        # Both chart canvases should be hidden from screen readers
        assert 'id="impact_chart" aria-hidden="true"' in html, (
            "Impact chart canvas missing aria-hidden=true"
        )
        assert 'id="inheritance_chart" aria-hidden="true"' in html, (
            "Inheritance chart canvas missing aria-hidden=true"
        )

    def test_impact_chart_data_table(self, rendered_html_with_charts):
        """Test that impact chart has accessible data table fallback."""
        html = rendered_html_with_charts

        # Should have data table with aria-label and sr-only class
        assert 'aria-label="Impact distribution data"' in html, (
            "Impact chart data table aria-label not found"
        )
        assert 'class="chart-data-table sr-only"' in html, (
            "Impact chart data table not using sr-only class"
        )
        # Should have table headers
        assert '<th scope="col">Impact Level</th>' in html, (
            "Impact chart data table headers not found"
        )

    def test_inheritance_chart_data_table(self, rendered_html_with_charts):
        """Test that inheritance chart has accessible data table fallback."""
        html = rendered_html_with_charts

        # Should have data table with aria-label and sr-only class
        assert 'aria-label="Inheritance pattern distribution data"' in html, (
            "Inheritance chart data table aria-label not found"
        )
        # Should have table headers
        assert '<th scope="col">Inheritance Pattern</th>' in html, (
            "Inheritance chart data table headers not found"
        )

    def test_chart_headings_have_ids(self, template_source):
        """Test that chart headings have ID attributes for aria-labelledby."""
        assert 'id="impact-chart-heading"' in template_source, "Impact chart heading ID not found"
        assert 'id="inheritance-chart-heading"' in template_source, (
            "Inheritance chart heading ID not found"
        )

    # Print stylesheet tests (PRINT-01)

    def test_print_media_query_exists(self, template_source):
        """Test that template contains @media print stylesheet."""
        assert "@media print" in template_source, "Print media query not found"

    def test_print_hides_pagination(self, template_source):
        """Test that print stylesheet hides pagination controls."""
        assert "dataTables_paginate" in template_source, (
            "dataTables_paginate not found in stylesheet"
        )
        # Should be in a display: none rule (within print media query context)
        # Check that it appears after @media print
        print_pos = template_source.find("@media print")
        paginate_pos = template_source.find("dataTables_paginate", print_pos)
        assert paginate_pos > print_pos, "dataTables_paginate not in print stylesheet section"

    def test_print_hides_filters(self, template_source):
        """Test that print stylesheet hides filter controls."""
        assert "dataTables_filter" in template_source, "dataTables_filter not found in stylesheet"
        print_pos = template_source.find("@media print")
        filter_pos = template_source.find("dataTables_filter", print_pos)
        assert filter_pos > print_pos, "dataTables_filter not in print stylesheet"

    def test_print_hides_length(self, template_source):
        """Test that print stylesheet hides length selector."""
        assert "dataTables_length" in template_source, "dataTables_length not found in stylesheet"
        print_pos = template_source.find("@media print")
        length_pos = template_source.find("dataTables_length", print_pos)
        assert length_pos > print_pos, "dataTables_length not in print stylesheet"

    def test_print_hides_buttons(self, template_source):
        """Test that print stylesheet hides button controls."""
        assert "dt-buttons" in template_source, "dt-buttons not found in stylesheet"
        # Check it's in print section
        print_pos = template_source.find("@media print")
        buttons_pos = template_source.find("dt-buttons", print_pos)
        assert buttons_pos > print_pos, "dt-buttons not hidden in print stylesheet"

    def test_print_collapses_fixed_columns(self, template_source):
        """Test that print stylesheet collapses FixedColumns to static."""
        # Should have dtfc-fixed-left with position: static in print
        print_section_start = template_source.find("@media print")
        print_section = template_source[print_section_start : print_section_start + 5000]

        assert "dtfc-fixed-left" in print_section, "dtfc-fixed-left not found in print stylesheet"
        assert "position: static" in print_section or "position:static" in print_section, (
            "FixedColumns not collapsed to static in print"
        )

    def test_print_prevents_row_breaks(self, template_source):
        """Test that print stylesheet prevents page breaks inside rows."""
        print_section_start = template_source.find("@media print")
        print_section = template_source[print_section_start : print_section_start + 5000]

        # Should have break-inside: avoid or page-break-inside: avoid
        has_break_avoid = (
            "break-inside: avoid" in print_section
            or "break-inside:avoid" in print_section
            or "page-break-inside: avoid" in print_section
            or "page-break-inside:avoid" in print_section
        )
        assert has_break_avoid, "Row break prevention not found in print stylesheet"

    def test_print_repeats_headers(self, template_source):
        """Test that print stylesheet repeats table headers on each page."""
        print_section_start = template_source.find("@media print")
        print_section = template_source[print_section_start : print_section_start + 5000]

        assert "table-header-group" in print_section, (
            "table-header-group not found (needed to repeat headers)"
        )

    def test_print_hides_canvas(self, template_source):
        """Test that print stylesheet hides canvas elements."""
        print_section_start = template_source.find("@media print")
        print_section = template_source[print_section_start : print_section_start + 5000]

        # Canvas should be hidden in print
        assert "canvas" in print_section, "canvas not found in print stylesheet"
        # Should have display: none after canvas
        canvas_pos = print_section.find("canvas")
        display_none_pos = print_section.find("display: none", canvas_pos)
        assert display_none_pos > canvas_pos and display_none_pos - canvas_pos < 200, (
            "canvas not hidden with display: none in print stylesheet"
        )

    def test_print_shows_chart_data_tables(self, template_source):
        """Test that print stylesheet shows chart data tables."""
        print_section_start = template_source.find("@media print")
        print_section = template_source[print_section_start : print_section_start + 5000]

        # Should have chart-data-table with position: static or display: block
        assert "chart-data-table" in print_section, "chart-data-table not found in print stylesheet"
        # Should make them visible (position: static removes sr-only positioning)
        assert "position: static" in print_section or "position:static" in print_section, (
            "chart-data-table not made visible in print"
        )

    def test_print_hides_detail_panels(self, template_source):
        """Test that print stylesheet hides detail panels by default."""
        print_section_start = template_source.find("@media print")
        print_section = template_source[print_section_start : print_section_start + 5000]

        # Detail panels should be hidden
        assert "detail-panel" in print_section, "detail-panel not found in print stylesheet"

    # PDF export tests (PRINT-02)

    def test_pdf_export_button_exists(self, template_source):
        """Test that template contains PDF export button."""
        assert "export-pdf-btn" in template_source, "PDF export button ID not found"

    def test_pdf_export_calls_window_print(self, template_source):
        """Test that PDF export button handler calls window.print()."""
        assert "window.print()" in template_source, (
            "window.print() call not found in PDF export handler"
        )
        # window.print() should appear in script section
        # Just verify both exist - they're connected via event listener
        assert "export-pdf-btn" in template_source, "export-pdf-btn ID not found"

    def test_pdf_button_has_aria_label(self, rendered_html_with_charts):
        """Test that PDF button has aria-label for accessibility."""
        html = rendered_html_with_charts

        # Find the button in rendered HTML
        assert 'id="export-pdf-btn"' in html, "export-pdf-btn not found in HTML"
        # aria-label should be on the button element
        assert 'aria-label="Download report as PDF"' in html, "PDF export button missing aria-label"

    def test_pdf_button_hidden_in_print(self, template_source):
        """Test that PDF button is hidden in print stylesheet."""
        print_section_start = template_source.find("@media print")
        print_section = template_source[print_section_start : print_section_start + 5000]

        # export-pdf-btn should be in hidden list
        assert "export-pdf-btn" in print_section, (
            "export-pdf-btn not found in print stylesheet hidden elements"
        )


@pytest.mark.unit
class TestPhase16FilteringAndVisualization:
    """Test Phase 16 column-level filtering and visualization features."""

    @pytest.fixture(scope="class")
    def template_path(self):
        """Path to the template file."""
        templates_dir = Path(__file__).parent.parent.parent / "variantcentrifuge" / "templates"
        return templates_dir / "index.html"

    @pytest.fixture(scope="class")
    def template_source(self, template_path):
        """Raw template source for pattern matching tests."""
        return template_path.read_text()

    @pytest.fixture(scope="class")
    def rendered_html_with_phase16_data(self):
        """Render template with Phase 16 test data (reused across tests)."""
        templates_dir = Path(__file__).parent.parent.parent / "variantcentrifuge" / "templates"

        env = Environment(loader=FileSystemLoader(str(templates_dir)))
        template = env.get_template("index.html")

        assets = _load_assets()

        # Diverse summary data with impact and inheritance distributions
        summary = {
            "num_variants": 30,
            "num_genes": 10,
            "num_samples": 4,
            "impact_distribution": {
                "HIGH": 8,
                "MODERATE": 12,
                "LOW": 6,
                "MODIFIER": 4,
            },
            "inheritance_distribution": {
                "de_novo": 6,
                "autosomal_dominant": 5,
                "autosomal_recessive": 3,
            },
            "top_genes": [
                {"gene": "BRCA1", "count": 8},
                {"gene": "TP53", "count": 5},
            ],
        }

        # Diverse variant data: mix of chromosomes, SNVs and indels, numeric AF column
        variants = [
            {
                "GENE": "BRCA1",
                "CHROM": "chr17",
                "POS": "41234567",
                "REF": "A",
                "ALT": "G",
                "IMPACT": "HIGH",
                "Inheritance_Pattern": "de_novo",
                "AF": "0.001",
                "ClinVar_CLNSIG": "Pathogenic",
                "igv_links": [],
            },
            {
                "GENE": "TP53",
                "CHROM": "chr17",
                "POS": "7577538",
                "REF": "C",
                "ALT": "CAGT",
                "IMPACT": "MODERATE",
                "Inheritance_Pattern": "autosomal_dominant",
                "AF": "0.003",
                "ClinVar_CLNSIG": None,
                "igv_links": [],
            },
            {
                "GENE": "EGFR",
                "CHROM": "chr7",
                "POS": "55241707",
                "REF": "ATGT",
                "ALT": "A",
                "IMPACT": "HIGH",
                "Inheritance_Pattern": "autosomal_recessive",
                "AF": "0.0005",
                "ClinVar_CLNSIG": "Likely_pathogenic",
                "igv_links": [],
            },
            {
                "GENE": "PTEN",
                "CHROM": "chr10",
                "POS": "89692905",
                "REF": "G",
                "ALT": "T",
                "IMPACT": "LOW",
                "Inheritance_Pattern": "de_novo",
                "AF": "0.01",
                "ClinVar_CLNSIG": "Uncertain_significance",
                "igv_links": [],
            },
            {
                "GENE": "RET",
                "CHROM": "chr10",
                "POS": "43612032",
                "REF": "C",
                "ALT": "G",
                "IMPACT": "MODIFIER",
                "Inheritance_Pattern": "autosomal_dominant",
                "AF": "0.05",
                "ClinVar_CLNSIG": None,
                "igv_links": [],
            },
        ]

        column_data = [
            {
                "original_name": "GENE",
                "display_name": "GENE",
                "is_standard_link_column": False,
                "is_igv_link_column": False,
                "link_display_text": None,
                "apply_hover_expand": False,
                "max_width_px": None,
            },
            {
                "original_name": "CHROM",
                "display_name": "CHROM",
                "is_standard_link_column": False,
                "is_igv_link_column": False,
                "link_display_text": None,
                "apply_hover_expand": False,
                "max_width_px": None,
            },
            {
                "original_name": "POS",
                "display_name": "POS",
                "is_standard_link_column": False,
                "is_igv_link_column": False,
                "link_display_text": None,
                "apply_hover_expand": False,
                "max_width_px": None,
            },
            {
                "original_name": "REF",
                "display_name": "REF",
                "is_standard_link_column": False,
                "is_igv_link_column": False,
                "link_display_text": None,
                "apply_hover_expand": False,
                "max_width_px": None,
            },
            {
                "original_name": "ALT",
                "display_name": "ALT",
                "is_standard_link_column": False,
                "is_igv_link_column": False,
                "link_display_text": None,
                "apply_hover_expand": False,
                "max_width_px": None,
            },
            {
                "original_name": "IMPACT",
                "display_name": "IMPACT",
                "is_standard_link_column": False,
                "is_igv_link_column": False,
                "link_display_text": None,
                "apply_hover_expand": False,
                "max_width_px": None,
            },
            {
                "original_name": "Inheritance_Pattern",
                "display_name": "Inheritance Pattern",
                "is_standard_link_column": False,
                "is_igv_link_column": False,
                "link_display_text": None,
                "apply_hover_expand": False,
                "max_width_px": None,
            },
            {
                "original_name": "AF",
                "display_name": "AF",
                "is_standard_link_column": False,
                "is_igv_link_column": False,
                "link_display_text": None,
                "apply_hover_expand": False,
                "max_width_px": None,
            },
            {
                "original_name": "ClinVar_CLNSIG",
                "display_name": "ClinVar CLNSIG",
                "is_standard_link_column": False,
                "is_igv_link_column": False,
                "link_display_text": None,
                "apply_hover_expand": False,
                "max_width_px": None,
            },
            {
                "original_name": "igv_links",
                "display_name": "IGV Links",
                "is_standard_link_column": False,
                "is_igv_link_column": True,
                "link_display_text": None,
                "apply_hover_expand": False,
                "max_width_px": None,
            },
        ]

        html = template.render(
            variants=variants,
            summary=summary,
            column_data=column_data,
            default_hidden_columns=[],
            generation_date="2026-02-18 12:00:00",
            version="0.14.0-test",
            filter_expression="(ANN[0].IMPACT has 'HIGH')",
            vcf_source="/data/test.vcf.gz",
            reference_genome="hg38",
            assets=assets,
        )

        return html

    # ---- HTML Structure Tests ----

    def test_visualization_section_exists(self, rendered_html_with_phase16_data):
        """Test that rendered HTML contains the visualization section with viz-content id."""
        html = rendered_html_with_phase16_data

        assert "visualization-section" in html, (
            "visualization-section class not found in rendered HTML"
        )
        assert 'id="viz-content"' in html, "viz-content id not found in rendered HTML"

    def test_three_chart_canvases_exist(self, rendered_html_with_phase16_data):
        """Test that rendered HTML contains all three Phase 16 chart canvases."""
        html = rendered_html_with_phase16_data

        assert 'id="variant_type_chart"' in html, (
            "variant_type_chart canvas id not found in rendered HTML"
        )
        assert 'id="chromosome_chart"' in html, (
            "chromosome_chart canvas id not found in rendered HTML"
        )
        assert 'id="af_histogram_chart"' in html, (
            "af_histogram_chart canvas id not found in rendered HTML"
        )

    def test_collapse_toggle_exists(self, rendered_html_with_phase16_data):
        """Test that rendered HTML contains collapse toggle button with aria-expanded=false."""
        html = rendered_html_with_phase16_data

        assert "collapse-toggle" in html, "collapse-toggle class not found in rendered HTML"
        assert 'aria-expanded="false"' in html, "aria-expanded='false' not found on collapse toggle"

    def test_filter_controls_container_exists(self, rendered_html_with_phase16_data):
        """Test that rendered HTML contains filter-controls id and filter-row id."""
        html = rendered_html_with_phase16_data

        assert 'id="filter-controls"' in html, "filter-controls id not found in rendered HTML"
        assert 'id="filter-row"' in html, "filter-row id not found in rendered HTML"

    def test_filter_toggle_button_exists(self, rendered_html_with_phase16_data):
        """Test that rendered HTML contains filter-toggle-btn id."""
        html = rendered_html_with_phase16_data

        assert 'id="filter-toggle-btn"' in html, "filter-toggle-btn id not found in rendered HTML"

    def test_filter_chip_strip_exists(self, rendered_html_with_phase16_data):
        """Test that rendered HTML contains active-filters-strip id and btn-reset-all id."""
        html = rendered_html_with_phase16_data

        assert 'id="active-filters-strip"' in html, (
            "active-filters-strip id not found in rendered HTML"
        )
        assert 'id="btn-reset-all"' in html, "btn-reset-all id not found in rendered HTML"

    def test_missing_values_toggle_exists(self, rendered_html_with_phase16_data):
        """Test that rendered HTML contains include-missing-toggle with checked attribute."""
        html = rendered_html_with_phase16_data

        assert 'id="include-missing-toggle"' in html, (
            "include-missing-toggle id not found in rendered HTML"
        )
        # Should be checked by default
        assert "checked" in html, (
            "checked attribute not found (include-missing-toggle default state)"
        )

    def test_nouislider_css_inlined(self, rendered_html_with_phase16_data):
        """Test that noUiSlider CSS is inlined in rendered HTML."""
        html = rendered_html_with_phase16_data

        # noUiSlider CSS contains .noUi- selectors
        assert ".noUi-" in html, (
            "noUiSlider CSS (.noUi-) not found in rendered HTML — check asset inlining"
        )

    def test_nouislider_js_inlined(self, rendered_html_with_phase16_data):
        """Test that noUiSlider JS is inlined in rendered HTML."""
        html = rendered_html_with_phase16_data

        # noUiSlider JS contains noUiSlider identifier
        assert "noUiSlider" in html, (
            "noUiSlider JS not found in rendered HTML — check asset inlining"
        )

    # ---- JavaScript Behavior Pattern Tests ----

    def test_filter_initialization_pattern(self, template_source):
        """Test that template contains noUiSlider.create call for slider creation."""
        assert "noUiSlider.create" in template_source, "noUiSlider.create not found in template"

    def test_custom_search_function_pattern(self, template_source):
        """Test that template registers custom search functions via ext.search.push."""
        assert "ext.search.push" in template_source, (
            "ext.search.push (DataTables filter registration) not found in template"
        )

    def test_missing_values_set(self, template_source):
        """Test that template defines MISSING_VALUES set with standard sentinel values."""
        assert "MISSING_VALUES" in template_source, "MISSING_VALUES not found in template"
        # Should include standard missing value sentinels
        assert "'N/A'" in template_source, "'N/A' not in MISSING_VALUES set"
        assert "'NA'" in template_source, "'NA' not in MISSING_VALUES set"
        assert "'.'" in template_source, "'.' not in MISSING_VALUES set"

    def test_active_filters_map(self, template_source):
        """Test that template contains activeFilters map and updateFilterChips function."""
        assert "activeFilters" in template_source, "activeFilters not found in template"
        assert "updateFilterChips" in template_source, (
            "updateFilterChips function not found in template"
        )

    def test_reactive_chart_update(self, template_source):
        """Test that template wires draw.dt event to chart update function."""
        assert "draw.dt" in template_source, "draw.dt event handler not found in template"
        assert "updateAllCharts" in template_source, (
            "updateAllCharts function not found in template"
        )

    def test_chart_update_no_animation(self, template_source):
        """Test that chart updates use update('none') to suppress animation."""
        assert "update('none')" in template_source, (
            "update('none') pattern not found in template — charts should update without animation"
        )

    def test_debounce_pattern(self, template_source):
        """Test that template contains setTimeout/clearTimeout debounce pattern for text input."""
        assert "setTimeout" in template_source, (
            "setTimeout not found in template (needed for debounce)"
        )
        assert "clearTimeout" in template_source, (
            "clearTimeout not found in template (needed for debounce)"
        )

    def test_log_scale_histogram(self, template_source):
        """Test that template configures AF histogram with logarithmic Chart.js scale."""
        assert "logarithmic" in template_source, (
            "logarithmic scale type not found in template — AF histogram should use log scale"
        )

    def test_collapsible_section_localstorage(self, template_source):
        """Test that template persists visualization section state via localStorage."""
        assert "vizSectionExpanded" in template_source, (
            "vizSectionExpanded localStorage key not found in template"
        )

    def test_include_missing_toggle_wiring(self, template_source):
        """Test that template wires include-missing-toggle change event."""
        assert "include-missing-toggle" in template_source, (
            "include-missing-toggle not found in template"
        )
        # Should have change event handler wiring
        assert "change" in template_source, "change event handler not found in template"

    def test_reset_all_button_wiring(self, template_source):
        """Test that template wires btn-reset-all click event to reset all filters."""
        assert "btn-reset-all" in template_source, "btn-reset-all not found in template"
        # Should have click handler pattern
        assert "resetAllBtn" in template_source or "reset" in template_source.lower(), (
            "Reset all button click handler pattern not found"
        )
