# VariantCentrifuge — Issue Assessment & Prioritized Roadmap

**Date:** 2026-02-12
**Assessed by:** Senior code review (Claude Code)
**Repo:** scholl-lab/variantcentrifuge (v0.8.25)
**Branch:** main (clean, 240280a)

---

## 1. Current State Assessment

### 1.1 Test Health

| Metric | Value |
|--------|-------|
| Tests collected | 930 |
| Collection errors | 16 (all same root cause) |
| Unit tests (marked) passing | 17/17 |
| Full suite | Cannot complete — some tests hang or run very long |

**Root cause of all 16 collection errors:**
```
ModuleNotFoundError: No module named 'psutil'
```
`psutil` is declared in `setup.py` `install_requires` but missing from `conda/environment.yml`.
Every test file that imports from `variantcentrifuge.stages` or `variantcentrifuge.cli` triggers this
because `processing_stages.py` line 23 does `import psutil` at module level.

**Unregistered pytest markers** (produce warnings):
- `performance` — used in `test_vectorized_replacer.py`
- `gene_burden` — used in `test_gene_burden_integration.py`

### 1.2 CI/CD Gap

Only one workflow exists: `.github/workflows/docs.yml` (Sphinx build + GitHub Pages deploy).

**Missing entirely:**
- Test workflow (pytest on push/PR)
- Lint workflow (black, flake8, isort checks)
- Release/publish workflow (PyPI, GitHub Releases)
- Container build workflow (Docker/GHCR)

This is the single most impactful infrastructure gap — no automated validation of any code change.

### 1.3 Packaging & Reproducibility

| Item | Status | Issue |
|------|--------|-------|
| `setup.py` | Functional | `python_requires=">=3.7"` is very old; 3.7 is EOL |
| `conda/environment.yml` | Broken | Named `annotation` (misleading), missing `psutil`/`jinja2`/`intervaltree`/`numpy` |
| Dockerfile | Does not exist | Issue #53 |
| Bioconda recipe | Does not exist | Issue #54 |
| `pyproject.toml` | Incomplete | Only has `[build-system]` and `[tool.black]`, no project metadata |

### 1.4 Code Quality Observations

- **Black/isort/flake8** configured and enforced via pre-commit — good.
- **Two pipeline architectures** coexist (classic + stage-based). The classic is default, stage-based is opt-in via `--use-new-pipeline`. This adds maintenance burden but is manageable.
- **Inheritance analysis** has a known correctness bug (#33) — comp_het false positives with multi-transcript variants.
- **Stage files are very large:** `analysis_stages.py` is ~137KB, `processing_stages.py` ~115KB. These are hard to review and maintain but not blocking.

---

## 2. Open Issues Inventory

| # | Title | Category | Labels | Assignee |
|---|-------|----------|--------|----------|
| 33 | comp_het false positives with split-snpeff-lines | Bug/Correctness | enhancement | berntpopp |
| 35 | Tumor filter for sample genotype combinations | Feature | enhancement | berntpopp |
| 47 | Add annotation hg38 liftover | Feature | enhancement | berntpopp |
| 52 | Annotation for same amino acid position pathogenic | Feature | enhancement | berntpopp |
| 53 | Docker version | Infrastructure | enhancement | berntpopp |
| 54 | Conda/mamba preset installation | Infrastructure | enhancement | berntpopp |
| 55 | Get filter options from annotation header | Feature | enhancement | berntpopp |
| 58 | Performance testing framework | Testing | — | — |
| 59 | Fix 16 failing tests + real-world docs | Bug/Testing | — | — |
| 60 | Real-world test datasets | Testing | — | — |
| 61 | Comprehensive report generation validation | Testing | — | — |
| 62 | Performance optimization (parallelization/chunking) | Performance | — | — |

**Total: 12 open issues** (1 correctness bug, 4 features, 2 infrastructure, 4 testing/perf, 1 mixed)

---

## 3. Prioritization Framework

Issues are prioritized by:
1. **Clinical correctness** — wrong results in a clinical genomics tool can affect patient care
2. **Developer velocity** — CI/test fixes unblock everything else
3. **User adoption** — Docker/conda make installation reliable
4. **Feature value** — clinical features expand the user base
5. **Long-term quality** — performance/testing infrastructure

---

## 4. Prioritized Roadmap

### Phase 1: Correctness & Test Foundation (1-2 days)

#### P1.1 — Fix #33: comp_het false positives with split-snpeff-lines

**Priority:** CRITICAL — analytical correctness bug
**Impact:** Can produce false positive compound heterozygous calls in clinical reports

**Problem:**
When `--split-snpeff-lines` produces multiple rows per variant (one per transcript),
and `--genotype-filter comp_het` is used without transcript filtering, a single het
variant appears as multiple rows. The comp_het logic counts rows, not unique variants,
producing false compound het calls.

**Solution:** Deduplicate by `(CHROM, POS, REF, ALT)` before het variant counting.

**Files to modify:**
- `variantcentrifuge/inheritance/comp_het.py` — add dedup before counting
- `variantcentrifuge/inheritance/comp_het_vectorized.py` — same fix in vectorized path
- `tests/test_inheritance/test_comp_het_split_lines.py` — new regression test

**Implementation sketch:**
```python
# In comp_het counting logic, before len(het_variants) >= 2 check:
unique_variants = het_variants_df.drop_duplicates(subset=["CHROM", "POS", "REF", "ALT"])
if len(unique_variants) >= 2:
    # proceed with compound het logic
```

**Validation:**
- Regression test: single variant with 3 transcripts -> must NOT be called comp_het
- Regression test: two variants with 3 transcripts each -> must be called comp_het
- Existing comp_het tests must still pass

**References:**
- [Best practices for variant calling in clinical sequencing (Genome Medicine, 2020)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7586657/)
- [Effective variant filtering in studies of rare human disease (PMC)](https://pmc.ncbi.nlm.nih.gov/articles/PMC8282602/)

---

#### P1.2 — Fix #59: Resolve test collection errors

**Priority:** HIGH — blocks all CI and developer confidence
**Actual scope:** Much smaller than the issue describes

**Root cause:** `psutil` not in conda env / dev install. All 16 errors are import failures.

**Tasks:**
1. Add missing dependencies to `conda/environment.yml`:
   - `psutil`
   - `jinja2`
   - `intervaltree`
   - `numpy`
2. Register missing pytest markers in `pytest.ini`:
   - `performance`
   - `gene_burden`
3. Verify all 930 tests collect and pass after `pip install -e .`
4. Investigate any tests that hang or take >60s — add `@pytest.mark.slow` if needed

**Files to modify:**
- `conda/environment.yml`
- `pytest.ini`

---

### Phase 2: CI/CD & Infrastructure (2-3 days)

#### P2.1 — Add CI test + lint workflow (NEW — not a filed issue)

**Priority:** HIGH — most impactful infrastructure gap
**Justification:** [AMP/CAP guidelines](https://www.sciencedirect.com/science/article/pii/S1525157817303732) require automated testing for clinical-grade bioinformatics pipelines.

**Implementation:**
```yaml
# .github/workflows/test.yml
name: Tests
on: [push, pull_request]
jobs:
  test:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: ["3.10", "3.12"]
    steps:
      - uses: actions/checkout@v4
      - uses: actions/setup-python@v5
        with:
          python-version: ${{ matrix.python-version }}
      - run: pip install -e ".[dev]"
      - run: black --check .
      - run: isort --check .
      - run: flake8 .
      - run: pytest -m "not slow" --tb=short
```

**Prerequisite:** Add `extras_require` to `setup.py`:
```python
extras_require={
    "dev": ["pytest", "pytest-cov", "pytest-mock", "black", "flake8",
            "flake8-docstrings", "isort", "pre-commit"],
}
```

---

#### P2.2 — Fix #54: Conda/mamba installation

**Priority:** MEDIUM
**Tasks:**
- Rename env `annotation` -> `variantcentrifuge`
- Add all missing Python deps (`psutil`, `jinja2`, `intervaltree`, `numpy`)
- Add `pip: - -e .` or `- variantcentrifuge` to install the package itself
- Test with `mamba env create -f conda/environment.yml && mamba activate variantcentrifuge && variantcentrifuge --help`

**Files to modify:**
- `conda/environment.yml`

---

#### P2.3 — Fix #53: Docker version

**Priority:** MEDIUM
**Justification:** [Containers in clinical bioinformatics](https://www.sciencedirect.com/science/article/pii/S1525157822000381) are essential for reproducibility and regulatory compliance.

**Implementation:**
```dockerfile
FROM continuumio/miniconda3:latest
COPY conda/environment.yml /tmp/environment.yml
RUN conda env create -f /tmp/environment.yml && conda clean -a
COPY . /app
WORKDIR /app
RUN conda run -n variantcentrifuge pip install -e .
ENTRYPOINT ["conda", "run", "--no-capture-output", "-n", "variantcentrifuge", "variantcentrifuge"]
```

**Additional tasks:**
- Add `.dockerignore` (exclude `.git`, `output/`, `__pycache__/`, `*.egg-info`)
- Add GitHub Actions workflow for container build + push to GHCR
- Document usage in README

---

### Phase 3: Clinical Feature Enhancements (1-2 weeks, parallelizable)

#### P3.1 — Fix #55: Get filter options from annotation header

**Priority:** MEDIUM — usability improvement
**Approach:**
- Add `--list-fields` CLI flag
- Parse VCF header via `bcftools view -h` or pysam
- Print available INFO/FORMAT fields with types and descriptions
- Optionally validate `--fields-to-extract` against available fields

**Files to modify:**
- `variantcentrifuge/cli.py` — add argument
- `variantcentrifuge/extractor.py` or new `variantcentrifuge/header_parser.py`

**Estimated scope:** ~150 lines

---

#### P3.2 — Fix #52: Same amino acid position pathogenic annotation

**Priority:** MEDIUM — high clinical value (supports ACMG PM5 criterion)
**Approach:**
- Parse `ANN[0].HGVS_P` to extract amino acid position
- Query ClinVar VCF (preloaded or on-the-fly) for pathogenic variants at same position
- Add annotation column `same_aa_pathogenic` (boolean or variant list)
- Implement as new annotation in `annotator.py` or as a new stage

**Dependencies:** Requires ClinVar VCF file path as input parameter

**Estimated scope:** ~200-300 lines

---

#### P3.3 — Fix #35: Tumor filter for sample genotype combinations

**Priority:** MEDIUM — extends oncology use case
**Approach:**
- Add `--tumor-normal-filter` option with preset patterns:
  - `somatic`: tumor `0/1` or `1/1`, normal `0/0`
  - `LOH`: tumor `1/1`, normal `0/1`
  - `germline`: both non-ref
- Require `--tumor-sample` and `--normal-sample` arguments
- Integrate as genotype filter preset alongside `comp_het`

**Files to modify:**
- `variantcentrifuge/filters.py` — new filter function
- `variantcentrifuge/cli.py` — new arguments
- `variantcentrifuge/config.json` — new presets

**Estimated scope:** ~200 lines

---

#### P3.4 — Fix #47: hg38 liftover annotation

**Priority:** MEDIUM — cross-reference capability
**Approach:**
- Optional integration with `CrossMap` or UCSC `liftOver`
- Add `--liftover-chain` argument
- Add `hg19_CHROM`, `hg19_POS` columns to output
- Cache chain file download

**External dependency:** Chain file (~300MB), CrossMap or liftOver binary

**Estimated scope:** ~200 lines

---

### Phase 4: Testing & Performance Infrastructure (ongoing)

#### P4.1 — Fix #58: Performance testing framework

**Priority:** LOW — enables validation of #62 but not blocking
**Approach:**
- Add `pytest-benchmark` to dev dependencies
- Benchmark critical paths: vectorized replacer, comp_het analysis, genotype replacement
- Establish baselines before any optimization work
- CI integration after basic CI exists (P2.1)

**Estimated scope:** ~500 lines, phased

---

#### P4.2 — Fix #62: Performance optimization

**Priority:** LOW — should wait for Phase 1 + Phase 2
**Justification:** Don't optimize code with known correctness bugs (#33). Don't optimize without benchmarks (#58) or CI (#P2.1) to validate changes.

**Approach:** Follow the phased roadmap in the issue:
- Phase 1: Dynamic chunk sizing, adaptive memory factors (after #58 baselines)
- Phase 2: Pipeline fusion, shared memory pools (after CI + benchmarks)
- Phase 3: Cache-aware layouts, SIMD (long-term)

---

#### P4.3 — Fix #60: Real-world test datasets

**Priority:** LOW — logistically complex
**Approach:**
- Expand existing `tests/fixtures/giab/` with GIAB reference samples
- Subset to manageable sizes (1000-5000 variants)
- Use `pooch` for lazy download with hash verification
- Avoid hosting PHI/clinical data

---

#### P4.4 — Fix #61: Comprehensive report validation

**Priority:** LOW — diminishing returns from full browser testing
**Approach:**
- HTML: BeautifulSoup structure validation (not Selenium)
- Excel: openpyxl sheet/column/row validation
- Snapshot testing for HTML output stability
- Skip full browser testing — fragile and slow

---

## 5. Dependency Graph

```
P1.1 (comp_het fix) ──────────────────────────────────┐
                                                       ├──> P4.2 (perf optimization)
P1.2 (test fixes) ──> P2.1 (CI workflow) ──> P4.1 ───┘
                       │
                       ├──> P2.2 (conda fix)
                       │
                       └──> P2.3 (Dockerfile) ──> needs P2.2

P3.1-P3.4 (features) can proceed independently in parallel

P4.3 (datasets) ──> P4.4 (report validation)
```

---

## 6. Quick Wins (can be done in < 1 hour each)

1. **Register missing pytest markers** in `pytest.ini` — 2 lines
2. **Add `psutil` to conda env** — 1 line
3. **Rename conda env** `annotation` -> `variantcentrifuge` — 1 line
4. **Add `.dockerignore`** — new file, ~10 lines
5. **Update `python_requires`** in `setup.py` from `>=3.7` to `>=3.10` (3.7-3.9 are EOL)

---

## 7. Risk Assessment

| Risk | Severity | Mitigation |
|------|----------|------------|
| comp_het false positives in clinical use (#33) | HIGH | Fix immediately, add regression tests |
| No CI — regressions go undetected | HIGH | Add test workflow before any feature work |
| Conda env missing deps — broken install | MEDIUM | Fix environment.yml, test in clean env |
| Performance work without benchmarks (#62) | MEDIUM | Defer until #58 baselines established |
| Large stage files (100KB+) hard to review | LOW | Accept for now; split only if adding features |
| Python 3.7 min version claim is misleading | LOW | Bump to 3.10, test in CI matrix |

---

## 8. References

- [AMP/CAP Standards for Validating NGS Bioinformatics Pipelines](https://www.sciencedirect.com/science/article/pii/S1525157817303732)
- [Best practices for variant calling in clinical sequencing](https://pmc.ncbi.nlm.nih.gov/articles/PMC7586657/)
- [Containers in Bioinformatics: Best Practices in Molecular Pathology](https://www.sciencedirect.com/science/article/pii/S1525157822000381)
- [Effective variant filtering in studies of rare human disease](https://pmc.ncbi.nlm.nih.gov/articles/PMC8282602/)
- [Bionitio: best practices for bioinformatics CLI software](https://academic.oup.com/gigascience/article/8/9/giz109/5572530)
- [CREDO: Reproducible Docker for bioinformatics](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-024-05695-9)
- [Benchmarking accelerated NGS analysis pipelines (2025)](https://academic.oup.com/bioinformaticsadvances/article/5/1/vbaf085/8132977)
