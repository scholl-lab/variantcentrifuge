# Agent Guidance Modernization Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use
> `superpowers:subagent-driven-development` (recommended) or
> `superpowers:executing-plans` to implement this plan task-by-task. Steps use
> checkbox (`- [ ]`) syntax for tracking.

**Goal:** Establish shared multi-LLM repository guidance and enforce the
600-line production-module budget locally and in CI.

**Architecture:** Root `AGENTS.md` becomes the concise, tool-neutral source of
truth. `CLAUDE.md` becomes a Claude-specific pointer. A standard-library LOC
checker reads a checked-in baseline allowlist, and Make, pre-commit, and GitHub
Actions invoke the same checker.

**Tech Stack:** Markdown, Python 3.10 standard library, Make, pre-commit,
GitHub Actions, pytest.

---

## Files

- Create: `AGENTS.md`
- Modify: `CLAUDE.md`
- Create: `.loc-allowlist`
- Create: `scripts/check_file_size.py`
- Create: `tests/unit/test_check_file_size.py`
- Modify: `Makefile`
- Modify: `.pre-commit-config.yaml`
- Modify: `.github/workflows/test.yml`

### Task 1: Add a tested LOC guard

- [ ] Add `tests/unit/test_check_file_size.py` with tests that create a
  600-line file, a 601-line file, and an allowlisted 601-line file in a temporary
  package. Assert the checker accepts the first and the matching allowlist case
  and rejects the unallowlisted oversized file.
- [ ] Run the focused tests and confirm they fail because the checker does not
  yet exist:

  ```bash
  python -m pytest tests/unit/test_check_file_size.py -q
  ```

- [ ] Implement `scripts/check_file_size.py` using only `argparse`, `pathlib`,
  and `sys`. It must recursively inspect Python files below a supplied source
  root, skip any path segment named `tests`, parse an optional
  `path:maximum_lines` allowlist, and report all violations before exiting 1.
- [ ] Generate `.loc-allowlist` from the current production-module baseline.
  Every entry must name a current file over 600 lines and cap it at its current
  physical line count.
- [ ] Re-run the focused test and the production guard:

  ```bash
  python -m pytest tests/unit/test_check_file_size.py -q
  python scripts/check_file_size.py variantcentrifuge --allowlist .loc-allowlist
  ```

### Task 2: Wire the guard into contributor and CI flows

- [ ] Add `lint-loc` to the Makefile. It must run:

  ```bash
  $(PYTHON) scripts/check_file_size.py variantcentrifuge --allowlist .loc-allowlist
  ```

- [ ] Include `lint-loc` in `ci-check` and update its progress count.
- [ ] Add a local pre-commit hook that invokes the same command with
  `pass_filenames: false`.
- [ ] Add a GitHub Actions `lint-loc` step in the existing lint job using
  `python scripts/check_file_size.py variantcentrifuge --allowlist .loc-allowlist`.
- [ ] Run `make lint-loc` and inspect the changed CI/pre-commit configuration.

### Task 3: Publish modern shared guidance

- [ ] Create `AGENTS.md` as the source of truth for Codex, Claude Code, Gemini,
  and human contributors. Include project layout, safe commands, validation,
  no-destructive-worktree rules, pipeline safety, test placement, dependencies,
  documentation, and the 600-line refactor-on-touch policy.
- [ ] Replace `CLAUDE.md` with a short pointer to `AGENTS.md`, retaining only
  Claude entrypoint details and explicit safety rules.
- [ ] Verify the guide references existing paths and Makefile targets.

### Task 4: Run the completion checks

- [ ] Run focused LOC tests:

  ```bash
  python -m pytest tests/unit/test_check_file_size.py -q
  ```

- [ ] Run the full local CI gate:

  ```bash
  make ci-check
  ```

- [ ] Inspect `git diff --check`, `git diff --stat`, and `git status --short`.
- [ ] Do not commit; this repository requires explicit user permission before
  committing.
