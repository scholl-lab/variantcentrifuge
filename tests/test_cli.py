# File: tests/test_cli.py
# Location: variantcentrifuge/tests/test_cli.py

"""
Tests for CLI module.

This file contains tests ensuring the CLI runs and shows help correctly.
"""

import json
import os
import subprocess
import tempfile
from unittest.mock import patch, MagicMock
import pytest


def test_cli_help():
    """Test that the CLI help message can be displayed."""
    cmd = ["python", "-m", "variantcentrifuge.cli", "--help"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    assert "usage:" in result.stdout


def test_bcftools_prefilter_in_help():
    """Test that the bcftools-prefilter argument is in the help text."""
    cmd = ["python", "-m", "variantcentrifuge.cli", "--help"]
    result = subprocess.run(cmd, capture_output=True, text=True)
    assert result.returncode == 0
    assert "--bcftools-prefilter" in result.stdout
    # Check that the help text mentions bcftools expression
    assert "bcftools expression" in result.stdout.lower()


class TestShowCheckpointStatus:
    """Test the --show-checkpoint-status functionality."""
    
    def test_show_checkpoint_status_no_state(self, tmp_path):
        """Test --show-checkpoint-status when no checkpoint state exists."""
        with patch('sys.argv', ['variantcentrifuge', '--show-checkpoint-status', '--output-dir', str(tmp_path)]):
            with patch('builtins.print') as mock_print:
                with pytest.raises(SystemExit) as exc_info:
                    from variantcentrifuge.cli import main
                    main()
                
                assert exc_info.value.code == 0
                # Should print that no checkpoint state was found
                mock_print.assert_any_call("Checking checkpoint status for original monolithic pipeline...")
                mock_print.assert_any_call(f"No checkpoint state found in {tmp_path}")
    
    def test_show_checkpoint_status_with_new_pipeline_flag(self, tmp_path):
        """Test --show-checkpoint-status with --use-new-pipeline flag."""
        with patch('sys.argv', ['variantcentrifuge', '--show-checkpoint-status', '--use-new-pipeline', '--output-dir', str(tmp_path)]):
            with patch('builtins.print') as mock_print:
                with pytest.raises(SystemExit) as exc_info:
                    from variantcentrifuge.cli import main
                    main()
                
                assert exc_info.value.code == 0
                # Should indicate it's checking the new pipeline
                mock_print.assert_any_call("Checking checkpoint status for new stage-based pipeline...")
                mock_print.assert_any_call(f"No checkpoint state found in {tmp_path}")
    
    def test_show_checkpoint_status_with_config_file(self, tmp_path):
        """Test --show-checkpoint-status with config file that specifies new pipeline."""
        # Create a config file with new pipeline enabled
        config_file = tmp_path / "config.json"
        config_data = {
            "use_new_pipeline_architecture": True,
            "reference": "GRCh38.99"
        }
        config_file.write_text(json.dumps(config_data))
        
        with patch('sys.argv', ['variantcentrifuge', '--show-checkpoint-status', '--config', str(config_file), '--output-dir', str(tmp_path)]):
            with patch('builtins.print') as mock_print:
                with pytest.raises(SystemExit) as exc_info:
                    from variantcentrifuge.cli import main
                    main()
                
                assert exc_info.value.code == 0
                # Should detect new pipeline from config
                mock_print.assert_any_call("Checking checkpoint status for new stage-based pipeline...")
                mock_print.assert_any_call(f"No checkpoint state found in {tmp_path}")
    
    def test_show_checkpoint_status_with_existing_state(self, tmp_path):
        """Test --show-checkpoint-status when checkpoint state exists."""
        from variantcentrifuge.checkpoint import PipelineState
        
        # Create a checkpoint state
        state = PipelineState(str(tmp_path))
        state.initialize({"test_config": "value"}, "1.0.0")
        state.start_step("test_step")
        state.complete_step("test_step")
        
        with patch('sys.argv', ['variantcentrifuge', '--show-checkpoint-status', '--output-dir', str(tmp_path)]):
            with patch('builtins.print') as mock_print:
                with pytest.raises(SystemExit) as exc_info:
                    from variantcentrifuge.cli import main
                    main()
                
                assert exc_info.value.code == 0
                # Should print the pipeline status summary
                mock_print.assert_any_call("Checking checkpoint status for original monolithic pipeline...")
                
                # Check that it printed checkpoint information
                print_calls = [call.args[0] for call in mock_print.call_args_list]
                summary_printed = any("Pipeline State Summary:" in call for call in print_calls)
                assert summary_printed
    
    def test_show_checkpoint_status_config_file_error(self, tmp_path):
        """Test --show-checkpoint-status with invalid config file."""
        # Create an invalid config file
        config_file = tmp_path / "invalid_config.json"
        config_file.write_text("invalid json content")
        
        with patch('sys.argv', ['variantcentrifuge', '--show-checkpoint-status', '--config', str(config_file), '--output-dir', str(tmp_path)]):
            with patch('builtins.print') as mock_print:
                with pytest.raises(SystemExit) as exc_info:
                    from variantcentrifuge.cli import main
                    main()
                
                assert exc_info.value.code == 0
                # Should fall back to original pipeline when config loading fails
                mock_print.assert_any_call("Checking checkpoint status for original monolithic pipeline...")
                mock_print.assert_any_call(f"No checkpoint state found in {tmp_path}")
    
    def test_show_checkpoint_status_flag_precedence(self, tmp_path):
        """Test that --use-new-pipeline flag takes precedence over config file."""
        # Create a config file with new pipeline disabled
        config_file = tmp_path / "config.json"
        config_data = {
            "use_new_pipeline_architecture": False,
            "reference": "GRCh38.99"
        }
        config_file.write_text(json.dumps(config_data))
        
        with patch('sys.argv', ['variantcentrifuge', '--show-checkpoint-status', '--use-new-pipeline', '--config', str(config_file), '--output-dir', str(tmp_path)]):
            with patch('builtins.print') as mock_print:
                with pytest.raises(SystemExit) as exc_info:
                    from variantcentrifuge.cli import main
                    main()
                
                assert exc_info.value.code == 0
                # Flag should override config file
                mock_print.assert_any_call("Checking checkpoint status for new stage-based pipeline...")
                mock_print.assert_any_call(f"No checkpoint state found in {tmp_path}")
    
    def test_show_checkpoint_status_logging_level(self, tmp_path):
        """Test --show-checkpoint-status with different logging levels."""
        with patch('sys.argv', ['variantcentrifuge', '--show-checkpoint-status', '--log-level', 'DEBUG', '--output-dir', str(tmp_path)]):
            with patch('builtins.print') as mock_print:
                with pytest.raises(SystemExit) as exc_info:
                    from variantcentrifuge.cli import main
                    main()
                
                assert exc_info.value.code == 0
                # Should work with different log levels
                mock_print.assert_any_call("Checking checkpoint status for original monolithic pipeline...")
                mock_print.assert_any_call(f"No checkpoint state found in {tmp_path}")
    
    def test_show_checkpoint_status_detailed_summary(self, tmp_path):
        """Test --show-checkpoint-status displays detailed checkpoint information."""
        from variantcentrifuge.checkpoint import PipelineState
        
        # Create a more complex checkpoint state
        state = PipelineState(str(tmp_path))
        state.initialize({"gene_name": "BRCA1", "vcf_file": "test.vcf"}, "1.0.0")
        
        # Add multiple steps
        state.start_step("step1")
        state.complete_step("step1")
        
        state.start_step("step2")
        state.complete_step("step2")
        
        state.start_step("step3")
        state.fail_step("step3", "Test error")
        
        with patch('sys.argv', ['variantcentrifuge', '--show-checkpoint-status', '--output-dir', str(tmp_path)]):
            with patch('builtins.print') as mock_print:
                with pytest.raises(SystemExit) as exc_info:
                    from variantcentrifuge.cli import main
                    main()
                
                assert exc_info.value.code == 0
                
                # Verify detailed information is printed
                print_calls = [call.args[0] for call in mock_print.call_args_list]
                full_output = "\n".join(print_calls)
                
                assert "Pipeline State Summary:" in full_output
                assert "step1" in full_output
                assert "step2" in full_output
                assert "step3" in full_output
                assert "✓" in full_output  # Completed steps
                assert "✗" in full_output  # Failed step
                assert "Test error" in full_output  # Error message
