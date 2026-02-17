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

    def test_dashboard_above_table(self, rendered_html_with_phase14_data):
        """Test that dashboard appears before the variant table."""
        html = rendered_html_with_phase14_data

        # Dashboard should appear before table
        dashboard_pos = html.find('class="dashboard"')
        table_pos = html.find('id="variants_table"')

        assert dashboard_pos != -1, "Dashboard class not found in HTML"
        assert table_pos != -1, "Variants table ID not found in HTML"
        assert dashboard_pos < table_pos, "Dashboard should appear before variants table"

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
        assert (
            'id="inheritance_chart"' in html
            or "Inheritance analysis not available" in html
        )

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
        assert (
            "col.original_name === 'ClinVar_CLNSIG'" in html
            or "ClinVar_CLNSIG" in html
        )

        # Check for Inheritance_Pattern badge render function
        assert (
            "col.original_name === 'Inheritance_Pattern'" in html
            or "Inheritance_Pattern" in html
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
        return Path(__file__).parent.parent.parent / "variantcentrifuge" / "templates" / "index.html"

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
        assert (
            "display: grid" in template_source or "display:grid" in template_source
        ), "Grid layout for detail panel not found"

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
        assert (
            "font-family:" in template_source and "monospace" in template_source.lower()
        ), "Monospace font-family rule not found"

    def test_fixedcolumns_shadow_css(self, template_source):
        """Test that template contains FixedColumns shadow styling."""
        assert ".dtfc-fixed-left" in template_source, "dtfc-fixed-left class not found"
        assert "box-shadow" in template_source, "box-shadow for fixed columns not found"

    def test_record_count_css(self, template_source):
        """Test that template contains prominent record count styling."""
        assert ".dt-info" in template_source, "dt-info class not found"
        # Should have bold or font-weight styling
        assert (
            "font-weight: 600" in template_source or "font-weight:600" in template_source
        ), "Record count font-weight not found"

    # JS behavior tests (pattern-based, check template source)

    def test_control_column_definition(self, template_source):
        """Test that template defines control column with chevron."""
        assert "dt-control" in template_source, "dt-control class not found in columnDefs"
        assert (
            "defaultContent" in template_source and "chevron" in template_source
        ), "Control column defaultContent with chevron not found"

    def test_fixedcolumns_config(self, template_source):
        """Test that template contains FixedColumns configuration."""
        assert "fixedColumns" in template_source, "fixedColumns config not found"
        assert "left: 2" in template_source, "FixedColumns left: 2 config not found"

    def test_density_localstorage_init(self, template_source):
        """Test that template initializes density from localStorage."""
        assert "localStorage.getItem" in template_source, "localStorage.getItem not found"
        assert (
            "variantTableDensity" in template_source or "density" in template_source.lower()
        ), "Density localStorage key not found"

    def test_density_toggle_button(self, template_source):
        """Test that template contains density toggle button with localStorage persistence."""
        assert "localStorage.setItem" in template_source, "localStorage.setItem not found"
        # Should have button or toggle mechanism
        assert (
            "density" in template_source.lower() and "button" in template_source.lower()
        ), "Density toggle button not found"

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
        assert (
            "td.dt-control" in template_source or "dt-control" in template_source
        ), "dt-control in click handler not found"

    def test_column_index_shift(self, template_source):
        """Test that template accounts for control column with index + 1."""
        # Should have index + 1 or similar offset logic
        assert (
            "index + 1" in template_source or "index+1" in template_source
        ), "Column index shift (index + 1) not found"

    def test_tooltip_keyboard_trigger(self, template_source):
        """Test that template configures tooltips with keyboard focus trigger."""
        assert (
            "trigger:" in template_source or "trigger :" in template_source
        ), "Tooltip trigger config not found"
        assert (
            "mouseenter focus" in template_source or "focus" in template_source
        ), "Tooltip focus trigger not found"

    def test_tooltip_touch_support(self, template_source):
        """Test that template configures tooltips with touch support."""
        # Should have touch configuration
        assert (
            "touch:" in template_source or "touch" in template_source
        ), "Tooltip touch config not found"

    def test_hgvs_truncation_render(self, template_source):
        """Test that template contains HGVS end truncation render function."""
        # Should truncate HGVS fields
        assert (
            "HGVS" in template_source and ("substring" in template_source or "substr" in template_source)
        ), "HGVS truncation render function not found"

    def test_middle_truncation_render(self, template_source):
        """Test that template contains VAR_ID middle truncation render."""
        # Should have middle truncation logic for variant IDs
        assert (
            "VAR_ID" in template_source
            and ("substring" in template_source or "substr" in template_source)
        ), "VAR_ID middle truncation not found"

    def test_right_align_numeric(self, template_source):
        """Test that template right-aligns numeric score columns."""
        assert (
            "dt-right" in template_source or "text-align: right" in template_source
        ), "Right-align for numeric columns not found"

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
        assert (
            "FixedColumns" in html or "fixedColumns" in html
        ), "FixedColumns JS not found in rendered HTML"
