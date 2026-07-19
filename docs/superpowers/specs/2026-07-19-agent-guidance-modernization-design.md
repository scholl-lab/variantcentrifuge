# Agent Guidance Modernization Design

## Goal

Make VariantCentrifuge straightforward and safe for contributors using Codex,
Claude Code, Gemini, and comparable tools. Establish one concise, authoritative
repository guide, retain a minimal Claude-specific entrypoint, and make the
600-line source-file policy executable locally and in GitHub Actions.

## Context

Sibling Python repositories use `AGENTS.md` as the tool-neutral source of
truth and keep `CLAUDE.md` as a short pointer. They pair clear validation
commands with an executable file-size policy. VariantCentrifuge currently has
a detailed Claude-only guide, an existing Makefile/CI workflow, and production
modules that predate a 600-line policy.

## Design

### Shared instruction source

Create root `AGENTS.md` as the authoritative guide for humans and all coding
agents. It will cover the project map, safe commands, validation expectations,
working rules, pipeline-specific safety boundaries, and documentation and
dependency conventions. It will explicitly state that `Makefile` targets and
`pyproject.toml` are the sources of truth for routine commands and dependency
metadata.

Replace `CLAUDE.md` with a compact Claude Code entrypoint that points to
`AGENTS.md`. It will retain only the no-direct-pipeline-execution warning and
the no-commit-without-user-approval rule, both of which remain relevant to
Claude-specific workflows.

### File-size policy

Production Python files in `variantcentrifuge/` and root Python entrypoints are
subject to a hard 600 physical-line budget. Tests, generated documentation,
and vendored or dependency content are outside this budget.

The repository has pre-existing modules over the limit. A checked-in
`.loc-allowlist` records each current oversized production path and its present
line count as a ceiling. The checker fails when a new non-allowlisted source
file exceeds 600 lines or an allowlisted file grows. A contributor who changes
an already oversized module must also decompose it into cohesive units before
adding behavior, keeping public facades stable where practical. Reducing an
allowlisted module removes or lowers its allowance.

### Enforcement

Add a dependency-free Python checker under `scripts/`, exposed through
`make lint-loc`. Add it to `make ci-check`, GitHub Actions, and local
pre-commit. This makes one command (`make ci-check`) representative of the
repository quality gates and catches violations before review.

## Non-goals

- Refactor the existing oversized modules in this documentation-focused
  change.
- Change VariantCentrifuge runtime behavior, public CLI options, dependencies,
  or release metadata.
- Add tool-specific instructions for every possible LLM client.

## Verification

- The checker succeeds against the current allowlisted baseline and fails for
  a synthetic over-limit file or an increased allowlist candidate.
- `make lint-loc` is part of `make ci-check`.
- GitHub Actions and pre-commit invoke the checker.
- `AGENTS.md` documents the policy and `CLAUDE.md` delegates to it.
