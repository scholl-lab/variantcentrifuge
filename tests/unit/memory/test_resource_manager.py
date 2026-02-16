"""
Unit tests for ResourceManager
"""

from unittest.mock import mock_open, patch

import pytest

from variantcentrifuge.memory import ResourceManager


@pytest.mark.unit
def test_memory_detection_cli_override(monkeypatch):
    """Test that CLI --max-memory-gb parameter takes precedence."""
    config = {"max_memory_gb": 16.0}
    rm = ResourceManager(config=config)
    assert rm.memory_gb == 16.0
    assert rm._get_memory_source() == "CLI"


@pytest.mark.unit
def test_memory_detection_slurm(monkeypatch):
    """Test SLURM_MEM_PER_NODE environment variable detection."""
    monkeypatch.setenv("SLURM_MEM_PER_NODE", "32768")  # 32GB in MB
    rm = ResourceManager()
    assert rm.memory_gb == 32.0
    assert rm._get_memory_source() == "SLURM"


@pytest.mark.unit
def test_memory_detection_pbs_gb(monkeypatch):
    """Test PBS memory detection with GB units."""
    monkeypatch.setenv("PBS_RESC_MEM", "16gb")
    rm = ResourceManager()
    assert rm.memory_gb == 16.0
    assert rm._get_memory_source() == "PBS"


@pytest.mark.unit
def test_memory_detection_pbs_mb(monkeypatch):
    """Test PBS memory detection with MB units."""
    monkeypatch.setenv("PBS_RESC_MEM", "16000mb")
    rm = ResourceManager()
    assert abs(rm.memory_gb - 15.625) < 0.01  # 16000/1024
    assert rm._get_memory_source() == "PBS"


@pytest.mark.unit
def test_memory_detection_cgroup_v1(monkeypatch, tmp_path):
    """Test cgroup v1 memory limit detection."""
    # Create fake cgroup v1 file
    cgroup_file = tmp_path / "memory.limit_in_bytes"
    cgroup_file.write_text("17179869184")  # 16GB in bytes

    with (
        patch("pathlib.Path.exists", return_value=True),
        patch("builtins.open", mock_open(read_data="17179869184")),
    ):
        rm = ResourceManager()
        # Should use cgroup detection
        assert rm._get_cgroup_memory_limit() == 16.0


@pytest.mark.unit
def test_memory_detection_cgroup_v2_max(monkeypatch, tmp_path):
    """Test cgroup v2 handles 'max' string (unlimited)."""
    with (
        patch("pathlib.Path.exists", return_value=True),
        patch("builtins.open", mock_open(read_data="max")),
    ):
        rm = ResourceManager()
        # Should return None for 'max' and fall back to psutil
        assert rm._get_cgroup_memory_limit() is None


@pytest.mark.unit
def test_memory_detection_cgroup_v2_numeric(monkeypatch, tmp_path):
    """Test cgroup v2 with numeric limit."""
    with (
        patch("pathlib.Path.exists", return_value=True),
        patch("builtins.open", mock_open(read_data="8589934592")),  # 8GB
    ):
        rm = ResourceManager()
        assert rm._get_cgroup_memory_limit() == 8.0


@pytest.mark.unit
def test_memory_detection_fallback(monkeypatch):
    """Test fallback to psutil when no environment variables set."""
    # Clear any environment variables
    monkeypatch.delenv("SLURM_MEM_PER_NODE", raising=False)
    monkeypatch.delenv("PBS_RESC_MEM", raising=False)

    rm = ResourceManager()
    # Should use psutil fallback, returns some positive value
    assert rm.memory_gb > 0
    assert rm._get_memory_source() == "psutil"


@pytest.mark.unit
def test_cpu_detection():
    """Test CPU detection returns positive integer."""
    rm = ResourceManager()
    assert isinstance(rm.cpu_cores, int)
    assert rm.cpu_cores > 0


@pytest.mark.unit
def test_auto_chunk_size_small_dataset():
    """Test chunk size calculation for small dataset."""
    rm = ResourceManager(config={"max_memory_gb": 16.0})
    # 100 variants, 10 samples - should fit entirely in memory
    chunk_size = rm.auto_chunk_size(total_items=100, num_samples=10)
    # Should return full dataset size for small data
    assert chunk_size == 100


@pytest.mark.unit
def test_auto_chunk_size_large_dataset():
    """Test chunk size calculation for large dataset."""
    rm = ResourceManager(config={"max_memory_gb": 8.0})
    # 1M variants, 5000 samples - requires chunking
    chunk_size = rm.auto_chunk_size(total_items=1_000_000, num_samples=5000)
    # Should return bounded chunk size
    assert 100 <= chunk_size <= 1_000_000
    assert chunk_size < 1_000_000  # Should chunk the data


@pytest.mark.unit
def test_auto_chunk_size_respects_bounds():
    """Test chunk size respects min/max bounds."""
    rm = ResourceManager(config={"max_memory_gb": 0.1})  # Very low memory
    # With very low memory, should still enforce minimum of 100 for calculated max_items
    chunk_size = rm.auto_chunk_size(total_items=10000, num_samples=1000)
    # Calculated max_items will be at least 100 (enforced by bounds)
    assert chunk_size >= 100


@pytest.mark.unit
def test_auto_workers_memory_constrained():
    """Test worker count when memory is the limiting factor."""
    rm = ResourceManager(config={"max_memory_gb": 4.0})
    # Each task needs 2GB, so max 2 workers with 4GB * 0.8 safety factor
    workers = rm.auto_workers(task_count=10, memory_per_task_gb=2.0)
    # Should be constrained by memory (4 * 0.8 / 2 = 1.6 -> 1 worker)
    assert workers >= 1
    assert workers <= 10


@pytest.mark.unit
def test_auto_workers_cpu_constrained():
    """Test worker count when CPU is the limiting factor."""
    rm = ResourceManager(config={"max_memory_gb": 64.0})
    # High memory, low memory per task
    workers = rm.auto_workers(task_count=100, memory_per_task_gb=0.5)
    # Should be constrained by CPU cores
    assert workers <= rm.cpu_cores
    assert workers >= 1


@pytest.mark.unit
def test_auto_workers_task_constrained():
    """Test worker count when task count is the limiting factor."""
    rm = ResourceManager(config={"max_memory_gb": 64.0})
    # Only 2 tasks total
    workers = rm.auto_workers(task_count=2, memory_per_task_gb=1.0)
    # Can't have more workers than tasks
    assert workers <= 2


@pytest.mark.unit
def test_should_parallelize_small():
    """Test parallelization decision for small dataset."""
    rm = ResourceManager(min_items_for_parallel=100)
    assert rm.should_parallelize(50) is False
    assert rm.should_parallelize(99) is False


@pytest.mark.unit
def test_should_parallelize_large():
    """Test parallelization decision for large dataset."""
    rm = ResourceManager(min_items_for_parallel=100)
    assert rm.should_parallelize(100) is True
    assert rm.should_parallelize(500) is True
    assert rm.should_parallelize(10000) is True


@pytest.mark.unit
def test_get_summary():
    """Test get_summary returns expected keys."""
    rm = ResourceManager(config={"max_memory_gb": 16.0})
    summary = rm.get_summary()

    # Check all expected keys present
    assert "memory_gb" in summary
    assert "memory_source" in summary
    assert "cpu_cores" in summary
    assert "cpu_source" in summary
    assert "memory_safety_factor" in summary
    assert "overhead_factor" in summary
    assert "min_items_for_parallel" in summary
    assert "safe_memory_gb" in summary

    # Check values are reasonable
    assert summary["memory_gb"] == 16.0
    assert summary["memory_source"] == "CLI"
    assert summary["cpu_cores"] > 0
    assert summary["memory_safety_factor"] == 0.80
    assert summary["overhead_factor"] == 3.0
    assert summary["safe_memory_gb"] == 16.0 * 0.80


@pytest.mark.unit
def test_custom_parameters():
    """Test ResourceManager with custom parameters."""
    rm = ResourceManager(
        config={"max_memory_gb": 32.0},
        overhead_factor=2.0,
        memory_safety_factor=0.9,
        min_items_for_parallel=500,
    )

    summary = rm.get_summary()
    assert summary["memory_gb"] == 32.0
    assert summary["overhead_factor"] == 2.0
    assert summary["memory_safety_factor"] == 0.9
    assert summary["min_items_for_parallel"] == 500
    assert summary["safe_memory_gb"] == 32.0 * 0.9

    # Test custom threshold
    assert rm.should_parallelize(400) is False
    assert rm.should_parallelize(500) is True


@pytest.mark.unit
def test_slurm_invalid_value(monkeypatch):
    """Test SLURM detection handles invalid values gracefully."""
    monkeypatch.setenv("SLURM_MEM_PER_NODE", "invalid")
    rm = ResourceManager()
    # Should fall back to psutil, not crash
    assert rm.memory_gb > 0


@pytest.mark.unit
def test_pbs_invalid_value(monkeypatch):
    """Test PBS detection handles invalid values gracefully."""
    monkeypatch.setenv("PBS_RESC_MEM", "invalid")
    rm = ResourceManager()
    # Should fall back to psutil, not crash
    assert rm.memory_gb > 0


# --- CPU detection tests ---


@pytest.mark.unit
def test_cpu_detection_slurm(monkeypatch):
    """Test SLURM_CPUS_PER_TASK environment variable detection."""
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "16")
    rm = ResourceManager()
    assert rm.cpu_cores == 16
    assert rm._get_cpu_source() == "SLURM"


@pytest.mark.unit
def test_cpu_detection_pbs(monkeypatch):
    """Test PBS_NUM_PPN environment variable detection."""
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)
    monkeypatch.setenv("PBS_NUM_PPN", "8")
    rm = ResourceManager()
    assert rm.cpu_cores == 8
    assert rm._get_cpu_source() == "PBS"


@pytest.mark.unit
def test_cpu_detection_slurm_invalid(monkeypatch):
    """Test SLURM CPU detection handles invalid values gracefully."""
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "invalid")
    rm = ResourceManager()
    # Should fall back to psutil/os, not crash
    assert rm.cpu_cores > 0


@pytest.mark.unit
def test_cpu_detection_pbs_invalid(monkeypatch):
    """Test PBS CPU detection handles invalid values gracefully."""
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)
    monkeypatch.setenv("PBS_NUM_PPN", "invalid")
    rm = ResourceManager()
    # Should fall back to psutil/os, not crash
    assert rm.cpu_cores > 0


@pytest.mark.unit
def test_cpu_detection_slurm_priority(monkeypatch):
    """Test SLURM takes priority over PBS for CPU detection."""
    monkeypatch.setenv("SLURM_CPUS_PER_TASK", "12")
    monkeypatch.setenv("PBS_NUM_PPN", "8")
    rm = ResourceManager()
    assert rm.cpu_cores == 12
    assert rm._get_cpu_source() == "SLURM"


@pytest.mark.unit
def test_cpu_detection_fallback(monkeypatch):
    """Test fallback to psutil/os when no environment variables set."""
    monkeypatch.delenv("SLURM_CPUS_PER_TASK", raising=False)
    monkeypatch.delenv("PBS_NUM_PPN", raising=False)
    rm = ResourceManager()
    assert rm.cpu_cores > 0
    assert rm._get_cpu_source() in ("psutil (physical)", "os (logical)")
