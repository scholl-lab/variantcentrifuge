# AGENTS.md — VariantCentrifuge

Shared instructions for human contributors and coding agents (Codex, Claude
Code, Gemini, and similar tools). This is the repository-wide source of truth;
keep client-specific files such as `CLAUDE.md` short and pointing here.

## Project

VariantCentrifuge is a research-use command-line pipeline for filtering,
extracting, and prioritising genetic variants from VCF files. It combines
bcftools/SnpEff/SnpSift processing, pandas analysis, inheritance analysis,
scoring, report generation, and optional Snakemake/HPC execution.

Treat genomic data as sensitive. Do not claim clinical validity or use the tool
as clinical decision support.

## Repository map

```text
variantcentrifuge/        application package
  pipeline_core/          stage framework, runner, context, workspace
  stages/                 setup, processing, analysis, and output stages
  inheritance/            inheritance and compound-heterozygous analysis
  association/            burden and association-testing components
tests/                    unit, integration, inheritance, performance, fixtures, mocks
config/ scoring/          filter and scoring configuration
workflow/ profiles/       Snakemake workflow and scheduler profiles
scripts/                  checked-in development and validation helpers
docs/source/              user and developer documentation source
```

Read the nearest relevant module and tests before changing behavior. Preserve
the stage dependency contract and pipeline artifact/workspace conventions.

## Commands

Prefer Make targets over ad-hoc command sequences:

```bash
make lint          # Ruff lint
make lint-loc      # enforce the production Python 600-line budget
make format-check  # Ruff formatting check
make typecheck     # mypy
make test-fast     # non-slow, non-integration tests
make ci-check      # local lint, LOC, format, type, and fast-test gate
make format        # apply Ruff formatting and safe lint fixes
```

Use `make install-dev` for a development install. In a fresh git worktree,
run `make build-ext` (the `make test*` targets do it for you): the Davies C
extension is untracked, and without it SKAT/COAST p-values silently fall back
to Liu and golden tests fail. Add or change
dependencies in `pyproject.toml`; do not hand-edit generated environment state.

## Safety and working rules

- **Never run `variantcentrifuge` directly** unless the user explicitly
  authorises the exact command and inputs. It can process sensitive data and
  run for hours. Prefer focused tests and mocked external tools.
- Do not delete, reset, revert, or overwrite changes you did not create.
- Keep changes task-scoped. Do not mix unrelated cleanup into a feature or fix.
- Do not commit, push, tag, or open pull requests without explicit user
  approval.
- Inspect `git status` and the relevant tests before editing. Use an isolated
  worktree for multi-file feature work when practical.
- Add regression and edge-case coverage for changed behavior, and improve test
  coverage on every touched module toward the **>90% coverage target**. Keep
  unit tests deterministic; integration tests must isolate or mock external
  bioinformatics tools.
- Keep new tests under the existing `tests/` structure; do not create a second
  test root.

## Engineering conventions

- Python 3.10+; use modern built-in generics (`list[str]`, `str | None`) and
  Ruff-compatible formatting.
- Keep configuration-driven behavior in the established JSON/configuration
  surfaces rather than hard-coding new filter or scoring rules.
- `variantcentrifuge/version.py` owns `__version__`; other code imports it.
- Retain the separation between setup, processing, analysis, and output stages.
  Add dependencies deliberately and keep stage inputs/outputs explicit.
- External tools (`bcftools`, `snpEff`, `SnpSift`, `bedtools`, and optionally
  `igv-reports`) are environmental dependencies. Test their integration through
  the repository mocks or controlled fixtures.
- Keep documentation source in `docs/source/`; do not edit generated
  `docs/build/` output.

## Code quality and coverage policy (>90% target)

Every change must improve both the code quality and test coverage of the modules
it touches:

- **Strive for >90% code coverage:** Aim for `>90%` coverage across the
  repository and enforce `>90%` coverage on all new modules, refactored splits,
  and modified functions/branches.
- **Ratchet quality on touched code:** When modifying a module, resolve nearby
  type gaps (`mypy --disallow-untyped-defs`), dead code, silent exception
  swallowing, or untested error/edge paths within the touched responsibility.
- **Test real behavior, not mocks:** Write deterministic unit tests covering
  normal execution, boundary inputs, and failure modes. Inspect uncovered lines
  with `pytest --cov=variantcentrifuge --cov-report=term-missing --cov-fail-under=0`.
- **Never lower the coverage floor:** `[tool.coverage.report] fail_under` in
  `pyproject.toml` is a minimum ratchet that may only move upward toward `90%+`.

## 600-line production-module policy

Production Python modules under `variantcentrifuge/` have a hard 600 physical
line budget. Tests are exempt.

- New production modules must be at or below 600 lines.
- `make lint-loc`, pre-commit, and GitHub Actions enforce the policy.
- Existing oversized modules are listed in `.loc-allowlist` at their current
  ceiling. They may shrink but must never grow.
- **When touching an allowlisted oversized module, refactor it into cohesive
  modules in the same change until each production module is at or below the
  budget.** Split by responsibility, not arbitrary line count; preserve public
  facades and add focused tests for the extracted behavior.
- If a module approaches 600 lines, design the split before adding behavior.

## Definition of done

Before claiming a code change is ready, run the narrowest relevant test with
coverage (`--cov=variantcentrifuge --cov-report=term-missing --cov-fail-under=0`) to verify the
touched module is thoroughly covered (striving for `>90%`), and then run
`make ci-check` when feasible. Inspect `git diff --check` and `git status
--short`. Report any checks you could not run and why.
