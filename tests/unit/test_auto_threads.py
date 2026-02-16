"""
Unit tests for --threads auto-detection and SnpSift memory flag.
"""

import pytest

from variantcentrifuge.filters import _snpsift_memory_flag


@pytest.mark.unit
def test_snpsift_memory_flag_with_known_memory():
    """Test SnpSift memory flag calculation with known memory."""
    cfg = {"max_memory_gb": 32.0}
    flag = _snpsift_memory_flag(cfg)
    # 32 * 0.25 = 8, clamped to [2, 16]
    assert flag == "-Xmx8g"


@pytest.mark.unit
def test_snpsift_memory_flag_low_memory():
    """Test SnpSift memory flag clamps to minimum 2GB."""
    cfg = {"max_memory_gb": 4.0}
    flag = _snpsift_memory_flag(cfg)
    # 4 * 0.25 = 1, clamped to min 2
    assert flag == "-Xmx2g"


@pytest.mark.unit
def test_snpsift_memory_flag_high_memory():
    """Test SnpSift memory flag clamps to maximum 16GB."""
    cfg = {"max_memory_gb": 128.0}
    flag = _snpsift_memory_flag(cfg)
    # 128 * 0.25 = 32, clamped to max 16
    assert flag == "-Xmx16g"


@pytest.mark.unit
def test_snpsift_memory_flag_no_config():
    """Test SnpSift memory flag with empty config uses auto-detection."""
    flag = _snpsift_memory_flag({})
    # Should return a valid -Xmx flag regardless of system
    assert flag.startswith("-Xmx")
    assert flag.endswith("g")
    # Extract numeric value
    gb = int(flag[4:-1])
    assert 2 <= gb <= 16


@pytest.mark.unit
def test_threads_auto_resolves_to_int(monkeypatch):
    """Test that --threads auto resolves to an integer CPU count."""
    from variantcentrifuge.memory.resource_manager import ResourceManager

    rm = ResourceManager()
    assert isinstance(rm.cpu_cores, int)
    assert rm.cpu_cores >= 1


@pytest.mark.unit
def test_threads_auto_respects_slurm(monkeypatch):
    """Test that --threads auto picks up SLURM_CPUS_PER_TASK."""
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "24")
    from variantcentrifuge.memory.resource_manager import ResourceManager

    rm = ResourceManager()
    assert rm.cpu_cores == 24
