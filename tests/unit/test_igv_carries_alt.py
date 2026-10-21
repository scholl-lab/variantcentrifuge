"""IGV reports must include split multi-allelic Mutect2 carriers (GT such as 0/1/0)."""

import pytest

from variantcentrifuge.generate_igv_report import _carries_alt


@pytest.mark.parametrize("gt", ["0/1", "1/1", "0|1", "1", "0/1/0", "0/1/0/0", "0|0|1"])
def test_carriers(gt):
    assert _carries_alt(gt)


@pytest.mark.parametrize("gt", ["0/0", "./.", ".", "0", "0/0/0", "0/./1", ""])
def test_non_carriers(gt):
    assert not _carries_alt(gt)
