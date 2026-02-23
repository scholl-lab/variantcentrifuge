# Phase 28: Tech Debt Cleanup - Research

**Researched:** 2026-02-22
**Domain:** Python class attributes, CLI argument wiring, ROADMAP documentation
**Confidence:** HIGH

## Summary

Phase 28 is a targeted cleanup of four specific tech debt items identified in the
v0.15.0-MILESTONE-AUDIT.md. The changes are all small and surgical: one missing class
attribute, three missing CLI arguments (plus one already-structured-for but not wired),
and one documentation typo fix. No new algorithms or abstractions are introduced.

All code patterns for the fixes already exist in the codebase. The new CLI args follow
the established pattern from `--skat-backend`, `--association-workers`, and similar flags.
The `parallel_safe` class attribute pattern is already implemented on every other test
class except RSKATTest.

**Primary recommendation:** Follow existing patterns exactly. Do not invent new patterns
for any of these items.

## Standard Stack

No new libraries required. All changes use existing infrastructure.

### Existing Patterns Being Replicated

| Pattern | Where to Copy From | Where to Apply |
|---------|-------------------|----------------|
| `parallel_safe: bool = False` | `allelic_series.py:277`, `r_backend.py:16` | `skat_r.py:RSKATTest` |
| CLI arg -> cfg[] -> _build_assoc_config | `cli.py:1156-1217` | Add 4 new cfg assignments |
| `add_argument()` in stats_group | `cli.py:423-428` (skat-backend pattern) | Add --skat-method, --min-cases, etc. |
| `_get(cli_key, json_key, default)` | `analysis_stages.py:2278` | Already exists for new fields |

## Architecture Patterns

### parallel_safe Attribute Pattern

All AssociationTest subclasses declare `parallel_safe` as a class-level typed attribute.
The engine uses `getattr(test, "parallel_safe", False)` — the default of False is correct
behavior but the docstring of skat_r.py (line 25) explicitly says the stage must declare it.

Current state of all test classes:
| Class | File | parallel_safe | Notes |
|-------|------|---------------|-------|
| FisherExactTest | fisher.py | `True` (declared) | Complete |
| LogisticBurdenTest | logistic_burden.py | `True` (declared) | Complete |
| LinearBurdenTest | linear_burden.py | `True` (declared) | Complete |
| PurePythonSKATTest | skat_python.py | `True` (declared) | Complete |
| COASTTest | allelic_series.py | `False` (declared) | Complete |
| PurePythonCOASTTest | allelic_series_python.py | `True` (declared) | Complete |
| **RSKATTest** | **skat_r.py** | **MISSING** | **Tech debt** |

The fix is a one-line addition at class level in RSKATTest (line ~82, after the class
docstring, before `def __init__`):
```python
parallel_safe: bool = False  # rpy2 restriction: main thread only
```

The docstring already says "The stage that uses RSKATTest must declare `parallel_safe=False`"
(skat_r.py line 25), which is slightly wrong — the test class itself should declare it, not
the stage. This is the exact fix.

Note: The ROADMAP success criterion 1 says "in both skat_r.py files" — there is only ONE
skat_r.py file at `variantcentrifuge/association/tests/skat_r.py`. The criterion wording
may be an artifact; the fix is just one file.

### CLI Argument -> Config Wiring Pattern

The established pattern for association CLI args:
1. Add `add_argument()` in `variantcentrifuge/cli.py` within the stats_group (around line 423)
2. Add `cfg["key"] = getattr(args, "key", default)` in the config assembly block (around line 1156)
3. The `_build_assoc_config_from_context()` in `analysis_stages.py` already has `_get()` calls
   for these fields; the CLI keys are already expected

#### skat_method Wiring

`_build_assoc_config_from_context` at line 2278:
```python
skat_method=_get("skat_method", default="SKAT", nullable=False),
```
The cli_key is `"skat_method"`. The CLI arg should be `--skat-method` and write
`cfg["skat_method"] = args.skat_method`.

Valid values: `"SKAT"`, `"Burden"`, `"SKATO"` (from `python_backend.py:692-699` and
`r_backend.py:450-483`). The python backend also accepts `"OPTIMAL.ADJ"` and `"OPTIMAL"`
as aliases for SKATO.

#### min_cases Wiring

`_build_assoc_config_from_context` at line 2280:
```python
min_cases=_get("association_min_cases", json_key="min_cases", default=200, nullable=False),
```
The cli_key is `"association_min_cases"`. CLI arg: `--min-cases` writing
`cfg["association_min_cases"] = args.min_cases`.

#### max_case_control_ratio Wiring

`_build_assoc_config_from_context` at lines 2281-2285:
```python
max_case_control_ratio=_get(
    "association_max_case_control_ratio",
    json_key="max_case_control_ratio",
    default=20.0,
    nullable=False,
),
```
The cli_key is `"association_max_case_control_ratio"`. CLI arg: `--max-case-control-ratio`
writing `cfg["association_max_case_control_ratio"] = args.max_case_control_ratio`.

#### min_case_carriers Wiring

`_build_assoc_config_from_context` at lines 2287-2290:
```python
min_case_carriers=_get(
    "association_min_case_carriers",
    json_key="min_case_carriers",
    default=10,
    nullable=False,
),
```
The cli_key is `"association_min_case_carriers"`. CLI arg: `--min-case-carriers` writing
`cfg["association_min_case_carriers"] = args.min_case_carriers`.

### ROADMAP Criterion Fix

The audit file `v0.15.0-MILESTONE-AUDIT.md` records (line 18):
```
"lambda_gc.tsv filename vs ROADMAP criterion specifying lambda_gc.txt (typo in criterion)"
```

**CURRENT STATE:** The ROADMAP Phase 22 criterion 3 (line 168) already reads:
```
3. The `--diagnostics-output` directory contains `lambda_gc.tsv` with one ...
```

The ROADMAP already says `.tsv`. However, the Phase 22 VERIFICATION.md (file
`.planning/phases/22-acat-o-and-diagnostics/22-VERIFICATION.md`) still records
the gap as unresolved. The Phase 28 success criterion 4 says "ROADMAP Phase 22
criterion correctly references `lambda_gc.tsv` (not `.txt`)".

**Planner conclusion:** ROADMAP is already correct. The fix required is to update the
Phase 22 VERIFICATION.md to mark the gap as resolved (ROADMAP now says .tsv). If the
planner prefers, criterion 4 can be verified by grep confirming `.tsv` in ROADMAP and
closing the verification gap documentation.

Implementation output (confirmed): `diagnostics.py:371` writes to `lambda_gc.tsv` (not
.txt). The implementation is correct. The ROADMAP criterion is already correct.

## Don't Hand-Roll

| Problem | Don't Build | Use Instead |
|---------|-------------|-------------|
| CLI arg parsing | Custom parsing | argparse `add_argument()` in existing group |
| Config precedence | New logic | `_get()` helper in `_build_assoc_config_from_context()` |
| JSON config validation | New validators | `VALID_ASSOCIATION_KEYS` frozenset + existing validator |

**Key insight:** All infrastructure for the new CLI args already exists in
`_build_assoc_config_from_context()`. The only gap is the CLI `add_argument()` calls and
the `cfg["key"] = args.key` assignments in cli.py.

## Common Pitfalls

### Pitfall 1: Wrong CLI Key Names for diagnostic thresholds

**What goes wrong:** Using `"min_cases"` as the cli_key in `_get()` instead of
`"association_min_cases"`.
**Why it happens:** The JSON config key (`min_cases`) differs from the cli_key
(`association_min_cases`) for the diagnostic thresholds. The `_get()` function takes
separate `cli_key` and `json_key` parameters for exactly this reason.
**How to avoid:** Use exactly the cli_keys already in `_build_assoc_config_from_context`:
- `"association_min_cases"` (not `"min_cases"`)
- `"association_max_case_control_ratio"` (not `"max_case_control_ratio"`)
- `"association_min_case_carriers"` (not `"min_case_carriers"`)

### Pitfall 2: skat_method Uses Plain Key (No association_ Prefix)

**What goes wrong:** Adding `cfg["association_skat_method"]` in cli.py instead of
`cfg["skat_method"]`.
**Why it happens:** The `_get()` call for skat_method at line 2278 uses cli_key
`"skat_method"` (no prefix), unlike the diagnostic thresholds which use the
`association_` prefix.
**How to avoid:** Match the exact cli_key in `_get()`. For skat_method:
`cfg["skat_method"] = getattr(args, "skat_method", "SKAT")`.

### Pitfall 3: Forgetting VALID_ASSOCIATION_KEYS for skat_method Validation

**What goes wrong:** Adding `--skat-method` CLI arg but not adding its validation to
`_validate_association_config_dict()`.
**Why it happens:** `skat_method` IS already in `VALID_ASSOCIATION_KEYS` (line 2052) and
in `str_keys` (line 2110). No change needed in the validator.
**How to avoid:** Verify before coding — check `VALID_ASSOCIATION_KEYS` at line 2044.
`skat_method` is already there.

### Pitfall 4: Missing Enum Validation for --skat-method

**What goes wrong:** Adding the CLI arg without constraining choices.
**How to avoid:** Use `choices=["SKAT", "Burden", "SKATO"]` in `add_argument()` OR use
`type=str` and validate in the existing `_validate_association_config_dict()`. Looking at
`--skat-backend`, it uses `choices=["auto", "r", "python"]` — follow this pattern for
`--skat-method`.

### Pitfall 5: parallel_safe in Wrong Location

**What goes wrong:** Adding `parallel_safe` as an instance attribute in `__init__` instead
of as a class attribute.
**Why it matters:** `getattr(test, "parallel_safe", False)` works on instances, so
instance attribute would also work. But all other tests declare it as a class attribute
for documentation clarity. Class attribute is the right pattern.
**How to avoid:** Place `parallel_safe: bool = False` at class body level, before
`__init__`, matching `allelic_series.py:277`.

### Pitfall 6: Tests for New CLI Args

**What goes wrong:** Writing tests that invoke cli.py directly and test arg parsing in
isolation.
**How to avoid:** Follow the pattern in `test_json_config.py` which uses
`_build_assoc_config_from_context()` with a mock context. For CLI arg propagation tests,
test that `cfg["skat_method"]` is set correctly after calling `args_to_config()` or test
`_build_assoc_config_from_context` with the new cli_keys set.

## Code Examples

### parallel_safe Class Attribute (verified pattern from allelic_series.py:277)

```python
class RSKATTest(AssociationTest):
    """..."""

    parallel_safe: bool = False  # rpy2 restriction: main thread only

    def __init__(self) -> None:
        ...
```

### CLI Argument Addition Pattern (verified from cli.py:423-428)

```python
stats_group.add_argument(
    "--skat-method",
    choices=["SKAT", "Burden", "SKATO"],
    default="SKAT",
    help=(
        "SKAT method variant: 'SKAT' (default), 'Burden' (burden-only), "
        "or 'SKATO' (SKAT-O omnibus). Requires SKAT test in --association-tests."
    ),
)
stats_group.add_argument(
    "--min-cases",
    type=int,
    default=None,
    help=(
        "Minimum number of cases threshold for diagnostics warning. "
        "Default: 200. Warn when n_cases < this value."
    ),
)
stats_group.add_argument(
    "--max-case-control-ratio",
    type=float,
    default=None,
    help=(
        "Maximum case:control ratio threshold for diagnostics warning. "
        "Default: 20.0. Warn when n_controls/n_cases > this value."
    ),
)
stats_group.add_argument(
    "--min-case-carriers",
    type=int,
    default=None,
    help=(
        "Minimum case carrier count per-gene threshold for diagnostics. "
        "Default: 10. Flag genes where case_carriers < this value."
    ),
)
```

### Config Assignment in cli.py (verified pattern from cli.py:1156)

```python
# skat_method — cli_key matches _get() call in _build_assoc_config_from_context
cfg["skat_method"] = getattr(args, "skat_method", "SKAT")

# Diagnostic thresholds — use association_ prefix to match _get() cli_keys
# None means "not set by CLI" so JSON can override (nullable=False fields need
# the key to be absent OR set to a value; for non-nullable, we write None only
# when the user didn't provide the flag)
_min_cases = getattr(args, "min_cases", None)
if _min_cases is not None:
    cfg["association_min_cases"] = _min_cases

_max_ratio = getattr(args, "max_case_control_ratio", None)
if _max_ratio is not None:
    cfg["association_max_case_control_ratio"] = _max_ratio

_min_carriers = getattr(args, "min_case_carriers", None)
if _min_carriers is not None:
    cfg["association_min_case_carriers"] = _min_carriers
```

**IMPORTANT NOTE on nullable behavior:** `_build_assoc_config_from_context` uses
`nullable=False` for these fields. The `_get()` helper for non-nullable fields checks
`if cli_key in cfg:` — so if the key is NOT in cfg, JSON wins. This means:
- For `skat_method`: always write `cfg["skat_method"]` (CLI always wins, safe to write
  default when not specified because the default matches `AssociationConfig` default)
- For diagnostic thresholds: only write cfg key when user explicitly passes the flag,
  so JSON config can still control it

However, looking at the existing pattern for `--pca-components` (always writes to cfg),
a simpler approach is to always write the value (using the argparse default). This is
what the existing code does for `skat_backend` (always writes, even with default "python").
Check which approach is preferred by looking at how `--association-workers` is handled:
`cfg["association_workers"] = getattr(args, "association_workers", 1)` — always written.

The safe choice is always write the cfg key (using argparse's `default=` value), matching
the existing pattern for all other association args.

## State of the Art

| Item | Current State | Target State |
|------|--------------|--------------|
| RSKATTest.parallel_safe | Missing (relies on getattr default) | `parallel_safe: bool = False` class attribute |
| --skat-method CLI arg | JSON-only via `association.skat_method` | CLI arg + JSON fallback |
| --min-cases CLI arg | JSON-only via `association.min_cases` | CLI arg + JSON fallback |
| --max-case-control-ratio | JSON-only via `association.max_case_control_ratio` | CLI arg + JSON fallback |
| --min-case-carriers | JSON-only via `association.min_case_carriers` | CLI arg + JSON fallback |
| ROADMAP lambda_gc criterion | Already says .tsv (line 168) | Verify + close verification gap |

## Open Questions

1. **skat_method default in argparse vs AssociationConfig**
   - What we know: `AssociationConfig.skat_method = "SKAT"` (default). `_build_assoc_config_from_context` uses `_get("skat_method", default="SKAT", nullable=False)`.
   - What's unclear: Should the argparse default be `"SKAT"` (always writes to cfg, CLI always wins) or `None` (only writes when user specifies, JSON can still override)?
   - Recommendation: Use `default="SKAT"` in argparse, always write cfg["skat_method"]. This matches how `skat_backend` is handled (`default="python"`, always written).

2. **ROADMAP lambda_gc typo: already fixed or still needs fix?**
   - What we know: ROADMAP line 168 currently says `.tsv`. The 22-VERIFICATION.md records it as a gap.
   - What's unclear: The Phase 28 criterion says "fix the typo" — but the ROADMAP already shows .tsv.
   - Recommendation: The planner should include a task to verify the ROADMAP says `.tsv` and update the 22-VERIFICATION.md gap to mark it resolved. The actual ROADMAP text requires no change.

3. **"both skat_r.py files" in success criterion 1**
   - What we know: There is only one `skat_r.py` in the codebase at `variantcentrifuge/association/tests/skat_r.py`.
   - What's unclear: The criterion says "in both skat_r.py files" — this may be a drafting error.
   - Recommendation: Add the attribute to the single existing `skat_r.py` file.

## Sources

### Primary (HIGH confidence)
- Direct file inspection of `variantcentrifuge/association/tests/skat_r.py` — RSKATTest class
- Direct file inspection of `variantcentrifuge/association/tests/skat_python.py` — PurePythonSKATTest
- Direct file inspection of `variantcentrifuge/association/tests/allelic_series.py` — COASTTest
- Direct file inspection of `variantcentrifuge/association/tests/allelic_series_python.py` — PurePythonCOASTTest
- Direct file inspection of `variantcentrifuge/association/tests/fisher.py` — FisherExactTest
- Direct file inspection of `variantcentrifuge/association/tests/logistic_burden.py` — LogisticBurdenTest
- Direct file inspection of `variantcentrifuge/association/tests/linear_burden.py` — LinearBurdenTest
- Direct file inspection of `variantcentrifuge/association/base.py` — AssociationConfig, AssociationTest
- Direct file inspection of `variantcentrifuge/association/engine.py` — getattr(test, "parallel_safe", False) at line 355
- Direct file inspection of `variantcentrifuge/cli.py` — association arg section, cfg assignments (lines 408-435, 1150-1217)
- Direct file inspection of `variantcentrifuge/stages/analysis_stages.py` — `_build_assoc_config_from_context`, `VALID_ASSOCIATION_KEYS` (lines 2044-2299)
- Direct file inspection of `.planning/ROADMAP.md` — Phase 22 criterion 3 (line 168), Phase 28 success criteria
- Direct file inspection of `.planning/v0.15.0-MILESTONE-AUDIT.md` — tech_debt section
- Direct file inspection of `.planning/phases/22-acat-o-and-diagnostics/22-VERIFICATION.md` — gap documentation

## Metadata

**Confidence breakdown:**
- RSKATTest parallel_safe fix: HIGH — single missing line, pattern fully established
- CLI arg wiring: HIGH — _build_assoc_config_from_context already has all _get() calls; only cli.py changes needed
- ROADMAP typo: HIGH — ROADMAP already says .tsv; fix is documentation/verification closure only
- Test strategy: HIGH — test_json_config.py shows exact pattern to follow

**Research date:** 2026-02-22
**Valid until:** Stable (no external dependencies; pure codebase knowledge)
