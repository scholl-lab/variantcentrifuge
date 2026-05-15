"""Tests for annotation config key normalization."""

from variantcentrifuge.config import normalize_annotation_config


def test_normalizes_canonical_bed_key_to_all_bed_aliases():
    cfg = {"annotate_bed_files": ["regions.bed"]}

    result = normalize_annotation_config(cfg)

    assert result is cfg
    assert result["annotate_bed_files"] == ["regions.bed"]
    assert result["annotate_bed"] == ["regions.bed"]


def test_normalizes_legacy_bed_key_to_canonical_key():
    cfg = {"annotate_bed": "regions.bed"}

    result = normalize_annotation_config(cfg)

    assert result["annotate_bed_files"] == ["regions.bed"]
    assert result["annotate_bed"] == ["regions.bed"]


def test_empty_legacy_bed_key_does_not_mask_canonical_bed_key():
    cfg = {"annotate_bed": [], "annotate_bed_files": ["regions.bed"]}

    result = normalize_annotation_config(cfg)

    assert result["annotate_bed_files"] == ["regions.bed"]
    assert result["annotate_bed"] == ["regions.bed"]


def test_normalizes_canonical_gene_list_key_to_all_gene_list_aliases():
    cfg = {"annotate_gene_lists": ["genes.txt"]}

    result = normalize_annotation_config(cfg)

    assert result["annotate_gene_lists"] == ["genes.txt"]
    assert result["annotate_gene_list_files"] == ["genes.txt"]
    assert result["annotate_gene_list"] == ["genes.txt"]


def test_normalizes_legacy_gene_list_key_to_canonical_key():
    cfg = {"annotate_gene_list": "genes.txt"}

    result = normalize_annotation_config(cfg)

    assert result["annotate_gene_lists"] == ["genes.txt"]
    assert result["annotate_gene_list_files"] == ["genes.txt"]
    assert result["annotate_gene_list"] == ["genes.txt"]


def test_normalizes_gene_list_files_key_to_canonical_key():
    cfg = {"annotate_gene_list_files": ["genes.txt"]}

    result = normalize_annotation_config(cfg)

    assert result["annotate_gene_lists"] == ["genes.txt"]
    assert result["annotate_gene_list_files"] == ["genes.txt"]
    assert result["annotate_gene_list"] == ["genes.txt"]


def test_empty_legacy_gene_list_key_does_not_mask_canonical_gene_list_key():
    cfg = {
        "annotate_gene_list": [],
        "annotate_gene_list_files": [],
        "annotate_gene_lists": ["genes.txt"],
    }

    result = normalize_annotation_config(cfg)

    assert result["annotate_gene_lists"] == ["genes.txt"]
    assert result["annotate_gene_list_files"] == ["genes.txt"]
    assert result["annotate_gene_list"] == ["genes.txt"]


def test_missing_annotation_keys_are_normalized_to_empty_lists():
    result = normalize_annotation_config({})

    assert result["annotate_bed_files"] == []
    assert result["annotate_bed"] == []
    assert result["annotate_gene_lists"] == []
    assert result["annotate_gene_list_files"] == []
    assert result["annotate_gene_list"] == []


def test_filters_empty_items_from_annotation_lists():
    cfg = {
        "annotate_bed_files": ["regions.bed", "", None],
        "annotate_gene_lists": ["genes.txt", "", None],
    }

    result = normalize_annotation_config(cfg)

    assert result["annotate_bed_files"] == ["regions.bed"]
    assert result["annotate_gene_lists"] == ["genes.txt"]


def test_treats_mapping_annotation_value_as_scalar_not_iterable_keys():
    cfg = {"annotate_bed_files": {"regions.bed": "repeat_region"}}

    result = normalize_annotation_config(cfg)

    assert result["annotate_bed_files"] == ["{'regions.bed': 'repeat_region'}"]
