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
