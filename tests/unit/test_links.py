"""Unit tests for build-aware external report links resolution and generation."""

import pytest

from variantcentrifuge.links import (
    DEFAULT_LINKS_BY_BUILD,
    add_links_to_table,
    normalize_reference,
    resolve_links,
)


@pytest.mark.unit
class TestNormalizeReference:
    """Tests for reference assembly normalization."""

    @pytest.mark.parametrize(
        "ref,expected",
        [
            ("GRCh37", "GRCh37"),
            ("GRCh37.75", "GRCh37"),
            ("GRCh37.p13", "GRCh37"),
            ("hg19", "GRCh37"),
            ("HG19", "GRCh37"),
            ("grch37", "GRCh37"),
            ("GRCh38", "GRCh38"),
            ("GRCh38.86", "GRCh38"),
            ("GRCh38.99", "GRCh38"),
            ("GRCh38.p13", "GRCh38"),
            ("hg38", "GRCh38"),
            ("HG38", "GRCh38"),
            ("grch38", "GRCh38"),
        ],
    )
    def test_supported_references_normalized(self, ref, expected):
        assert normalize_reference(ref) == expected

    @pytest.mark.parametrize("ref", ["mm10", "GRCm38", "canFam3", "unknown", "", None])
    def test_unsupported_references_return_none(self, ref):
        assert normalize_reference(ref) is None


@pytest.mark.unit
class TestResolveLinks:
    """Tests for link template resolution based on genome build."""

    def test_default_grch37_resolution(self):
        config = {"reference": "GRCh37.p13"}
        links = resolve_links(config)
        assert "Varsome" in links
        assert "hg19" in links["Varsome"]
        assert "dataset=gnomad_r2_1" in links.get("gnomAD_2", "")
        assert "hg=19" in links["SpliceAI"]
        assert config["links"] == links

    def test_grch38_resolution(self):
        config = {"reference": "GRCh38.99"}
        links = resolve_links(config)
        assert "Varsome" in links
        assert "hg38" in links["Varsome"]
        assert "dataset=gnomad_r4" in links.get("gnomAD_4", "")
        assert "hg=38" in links["SpliceAI"]
        assert "hg38" in links["Franklin"]
        assert "hg38" in links["autopvs1"]
        assert config["links"] == links

    def test_unsupported_build_produces_no_bundled_links(self):
        """Unsupported builds must NOT silently produce hg19 links."""
        config = {"reference": "mm10"}
        links = resolve_links(config)
        assert links == {}
        assert config["links"] == {}

    def test_explicit_custom_templates_preserved(self):
        custom = {"CustomPortal": "https://portal.example.org/{CHROM}-{POS}"}
        config = {"reference": "GRCh38", "links": custom}
        links = resolve_links(config)
        assert links == custom
        assert config["links"] == custom

    def test_grch38_overrides_stale_bundled_grch37_links(self):
        """When reference is GRCh38 but config contains stale default GRCh37 links."""
        config = {
            "reference": "GRCh38",
            "links": dict(DEFAULT_LINKS_BY_BUILD["GRCh37"]),
        }
        links = resolve_links(config)
        assert "hg38" in links["Varsome"]
        assert "dataset=gnomad_r4" in links["gnomAD_4"]


@pytest.mark.unit
class TestAddLinksToTableBuildAware:
    """Test URL generation with build-aware links."""

    def test_grch38_coordinates_generate_grch38_urls(self):
        table_lines = [
            "CHROM\tPOS\tREF\tALT\tID",
            "1\t55039909\tTC\tT\trs123",
        ]
        config = {"reference": "GRCh38.99"}
        link_configs = resolve_links(config)
        result = add_links_to_table(table_lines, link_configs)

        header = result[0].split("\t")
        row = result[1].split("\t")
        col_map = dict(zip(header, row, strict=False))

        assert col_map["Varsome"] == "https://varsome.com/variant/hg38/1-55039909-TC-T"
        assert "dataset=gnomad_r4" in col_map["gnomAD_4"]
        assert "&hg=38&" in col_map["SpliceAI"]
        assert col_map["Franklin"] == (
            "https://franklin.genoox.com/clinical-db/variant/snp/1-55039909-TC-T-hg38"
        )
        assert col_map["autopvs1"] == "https://autopvs1.bgi.com/variant/hg38/1-55039909-TC-T"
