"""Early-exit and interactive listing/resume helpers for VariantCentrifuge CLI."""

from __future__ import annotations

import argparse
import logging
import sys
from typing import Any

logger = logging.getLogger("variantcentrifuge")


def handle_listing_and_resume_options(args: argparse.Namespace, cfg: dict[str, Any]) -> None:
    """Handle early-exit listing options and interactive resume."""
    if getattr(args, "list_stages", False):
        from .display_utils import display_available_stages
        from .pipeline import create_stages_from_config
        from .stages.stage_registry import initialize_registry

        initialize_registry()
        try:
            stages = create_stages_from_config(cfg)
            display_available_stages(stages, cfg)
        except Exception as e:
            logger.error(f"Failed to create stages: {e}")
            sys.exit(1)
        sys.exit(0)

    if getattr(args, "list_checkpoints", False):
        from .checkpoint import PipelineState
        from .display_utils import display_enhanced_status

        pipeline_state = PipelineState(args.output_dir)
        if pipeline_state.load():
            display_enhanced_status(pipeline_state, [])
        else:
            print(f"❌ No checkpoint file found in {args.output_dir}")
            print("   Start a new pipeline run with --enable-checkpoint to create checkpoints.")
        sys.exit(0)

    if getattr(args, "interactive_resume", False):
        from .checkpoint import PipelineState
        from .interactive_resume import handle_interactive_resume
        from .pipeline import create_stages_from_config
        from .stages.stage_registry import initialize_registry

        initialize_registry()
        pipeline_state = PipelineState(args.output_dir)
        if not pipeline_state.load():
            print(f"❌ No checkpoint file found in {args.output_dir}")
            print("   Start a new pipeline run with --enable-checkpoint to create checkpoints.")
            sys.exit(1)

        try:
            stages = create_stages_from_config(cfg)
            stage_names = [stage.name for stage in stages]
            selected_stage = handle_interactive_resume(args, pipeline_state, stage_names)
            if selected_stage:
                cfg["resume_from"] = selected_stage
                logger.info(f"Interactive resume: Selected stage '{selected_stage}'")
            else:
                print("❌ No stage selected. Exiting.")
                sys.exit(0)
        except Exception as e:
            logger.error(f"Interactive resume failed: {e}")
            sys.exit(1)
