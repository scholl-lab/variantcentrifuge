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

    def test_grch38_per_entry_resolution_with_custom_link(self):
        """Default links upgrade to GRCh38 while preserving custom entries (#144)."""
        custom_url = "https://example.org/{CHROM}:{POS}"
        config = {
            "reference": "GRCh38.99",
            "links": {
                **DEFAULT_LINKS_BY_BUILD["GRCh37"],
                "CustomPortal": custom_url,
            },
        }
        links = resolve_links(config)
        assert "hg38" in links["Varsome"]
        assert "hg38" in links["Franklin"]
        assert "hg=38" in links["SpliceAI"]
        assert "dataset=gnomad_r4" in links["gnomAD_4"]
        assert "gnomAD_2" not in links
        assert links["CustomPortal"] == custom_url
        assert config["links"] == links

    def test_grch38_omits_bundled_links_when_absent_from_build_defaults(self):
        """When an unmodified bundled link is absent from target build defaults, omit it (#144)."""
        grch38_without_autopvs1_or_gnomad = {
            "ClinVar": DEFAULT_LINKS_BY_BUILD["GRCh38"]["ClinVar"],
            "Varsome": DEFAULT_LINKS_BY_BUILD["GRCh38"]["Varsome"],
            "Franklin": DEFAULT_LINKS_BY_BUILD["GRCh38"]["Franklin"],
            "SpliceAI": DEFAULT_LINKS_BY_BUILD["GRCh38"]["SpliceAI"],
        }
        config = {
            "reference": "GRCh38",
            "links_by_build": {"GRCh38": grch38_without_autopvs1_or_gnomad},
            "links": dict(DEFAULT_LINKS_BY_BUILD["GRCh37"]),
        }
        links = resolve_links(config)
        assert "autopvs1" not in links
        assert "gnomAD_2" not in links
        assert "gnomAD_4" not in links
        assert "Varsome" in links
        assert "hg38" in links["Varsome"]

    def test_custom_link_preserved_when_in_grch37_links_by_build_but_omitted_from_grch38(self):
        """Custom link in links and links_by_build['GRCh37'] is preserved on GRCh38 (#144)."""
        portal_url = "https://portal.example.org/{CHROM}-{POS}"
        config = {
            "reference": "GRCh38",
            "links": {
                **DEFAULT_LINKS_BY_BUILD["GRCh37"],
                "Portal": portal_url,
            },
            "links_by_build": {
                "GRCh37": {"Portal": portal_url},
                "GRCh38": {"Varsome": DEFAULT_LINKS_BY_BUILD["GRCh38"]["Varsome"]},
            },
        }
        links = resolve_links(config)
        assert "Portal" in links
        assert links["Portal"] == portal_url
        assert "Varsome" in links
        assert "hg38" in links["Varsome"]

    def test_custom_url_on_bundled_link_name_preserved_with_shipped_config(self):
        """Custom URL on bundled link name (Varsome) preserved with shipped config.json (#144)."""
        import json
        from pathlib import Path

        cfg_path = Path(__file__).parent.parent.parent / "variantcentrifuge" / "config.json"
        with open(cfg_path, encoding="utf-8") as f:
            config = json.load(f)

        config["reference"] = "GRCh38.99"
        custom_url = "https://custom.example.org/{CHROM}/{POS}"
        config["links"]["Varsome"] = custom_url

        resolved = resolve_links(config)
        assert resolved["Varsome"] == custom_url
        assert "hg=38" in resolved["SpliceAI"]
        assert "hg38" in resolved["Franklin"]
        assert "gnomad_r4" in resolved.get("gnomAD_4", "")

    def test_custom_links_by_build_not_overwritten_for_custom_portal(self):
        """Custom URL on custom link name (Portal) preserved even if in links_by_build (#144)."""
        config = {
            "reference": "GRCh38",
            "links": {"Portal": "https://grch37.custom/item"},
            "links_by_build": {
                "GRCh37": {"Portal": "https://grch37.custom/item"},
                "GRCh38": {"Portal": "https://grch38.custom/item"},
            },
        }
        resolved = resolve_links(config)
        assert resolved["Portal"] == "https://grch37.custom/item"

    @pytest.mark.parametrize(
        "links_dict",
        [
            {
                "gnomAD_4": "https://custom.gnomad.org/{CHROM}-{POS}",
                "gnomAD_2": DEFAULT_LINKS_BY_BUILD["GRCh37"]["gnomAD_2"],
            },
            {
                "gnomAD_2": DEFAULT_LINKS_BY_BUILD["GRCh37"]["gnomAD_2"],
                "gnomAD_4": "https://custom.gnomad.org/{CHROM}-{POS}",
            },
        ],
    )
    def test_custom_gnomad_preserved_regardless_of_dict_key_order(self, links_dict):
        """Explicit custom gnomAD_4 link is preserved and not overwritten by gnomAD_2 (#144)."""
        config = {
            "reference": "GRCh38",
            "links": dict(links_dict),
        }
        resolved = resolve_links(config)
        assert resolved["gnomAD_4"] == "https://custom.gnomad.org/{CHROM}-{POS}"

    def test_load_config_with_custom_destination_in_links_by_build_resolves_portal(self):
        """Custom destination added to links_by_build['GRCh38'] resolves for GRCh38 (#144)."""
        from variantcentrifuge.config import load_config

        cfg = load_config()
        cfg["reference"] = "GRCh38.99"
        portal_url = "https://portal.example/{CHROM}:{POS}:{REF}:{ALT}"
        cfg["links_by_build"]["GRCh38"]["Portal"] = portal_url

        resolved = resolve_links(cfg)
        assert "Portal" in resolved
        assert resolved["Portal"] == portal_url
        assert "Varsome" in resolved
        assert "hg38" in resolved["Varsome"]
        assert "gnomAD_4" in resolved


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
