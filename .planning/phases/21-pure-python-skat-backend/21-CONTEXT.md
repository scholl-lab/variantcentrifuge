# Phase 21: Pure Python SKAT Backend - Context

**Gathered:** 2026-02-20
**Status:** Ready for planning

<domain>
## Phase Boundary

Pure Python reimplementation of SKAT and SKAT-O that matches R output within tiered tolerance thresholds, using Davies ctypes (via CFFI) for exact p-values with Kuonen saddlepoint and Liu moment-matching fallbacks. Enables R-free operation. The R backend (Phase 20) serves as the correctness oracle for validation.

</domain>

<decisions>
## Implementation Decisions

### Fallback chain behavior
- Full three-tier chain: Davies -> Kuonen saddlepoint -> Liu moment-matching
- Matches modern best practice (SKATh, Wu et al. 2016; GENESIS; GMMAT)
- Each fallback logged at INFO level (always visible in logs, not just verbose mode)
- GENESIS+GMMAT hybrid triggering conditions for maximum safety:
  1. `ifault > 0` (Davies explicitly failed)
  2. `p >= 1.0` (mathematically invalid)
  3. `p < 1000 * machine_epsilon` (~2.2e-13, GENESIS-style underflow guard)
  4. `p <= 1e-5` (GMMAT-style proactive saddlepoint, catches false convergence zone)
- Davies parameters: `acc=1e-9`, `lim=100_000` (tight tolerance per SKATh recommendation)
- Saddlepoint implementation: Kuonen (1999) Lugannani-Rice formula, ~30 lines of math, Newton-Raphson for saddlepoint equation
- Liu as last resort only — known to be anti-conservative at genome-wide significance (25x type I error inflation at alpha=1e-8)

### Compilation & portability
- Build backend switch: hatchling -> setuptools (hatchling has no C extension support)
- CFFI API mode for qfc.c compilation (argon2-cffi / cryptography pattern)
  - `_build.py` declares C function signature, points to bundled qfc.c
  - CFFI builds CPython extension module at wheel-build time
  - `OptionalBuildExt` class makes compilation non-fatal from sdist
- qfc.c source: extracted from CompQuadForm CRAN package, standalone (no R dependencies)
- Compiled extension lives in package install directory (standard Python packaging)
- `try/except ImportError` at runtime falls back to Liu if extension unavailable
- `VARIANTCENTRIFUGE_NO_C_EXT=1` environment variable to skip build
- Target platforms: Linux gcc, macOS gcc/clang (Windows via WSL, no MinGW/MSVC)
- cffi added as runtime dependency (stable, widely installed in scientific Python)

### Validation tolerance & edge cases
- Tiered tolerance thresholds (tighter than roadmap's original 10%):
  - For p > 1e-4: require `|py_p - r_p| / r_p < 1e-4` (relative)
  - For p <= 1e-4: require `|log10(py_p) - log10(r_p)| < 0.05` (~12% relative on p-value)
- Three-phase validation strategy:
  1. Component-level: GLM residuals, eigenvalues, Davies p-value each compared against R via rpy2 with per-step tolerances
  2. Golden files: fixed test dataset with all R intermediate values saved as fixtures (no R needed in CI)
  3. Live rpy2 comparison: optional slow test that runs both backends on same data
- Critical implementation details for maximum R agreement:
  - `scipy.linalg.eigh(driver='evr')` to match R's DSYEVR
  - Extract response residuals from statsmodels GLM (not deviance/Pearson)
  - Eigenvalue filtering: `lambda > mean(lambda) / 100000` (match R exactly)
  - SKAT-O rho grid: `[0, 0.01, 0.04, 0.09, 0.25, 0.5, 1.0]` (match R exactly)
  - `scipy.integrate.quad(epsabs=1.49e-8)` for SKAT-O integration (matches R QUADPACK)
- Degenerate genes (rank < 2): include in output with p_value=NA + skip_reason column
- Zero-variant genes: same treatment (p_value=NA, skip_reason='zero_variants')

### Output & transparency
- `p_method` column: `"davies"` | `"saddlepoint"` | `"liu"` (mathematical names, per-gene)
- `p_converged` column: boolean — True if primary method (Davies) succeeded, False if fallback was needed
- `skip_reason` column: populated for genes with p_value=NA — values like `rank_deficient`, `zero_variants`, `insufficient_carriers`
- `skat_o_rho` column: optimal rho value for SKAT-O (0=pure SKAT, 1=pure burden), matches R backend output

### Claude's Discretion
- Exact saddlepoint Taylor expansion near q = E[Q] (singularity handling)
- Newton-Raphson convergence criteria for saddlepoint equation
- CFFI build script structure and module naming
- pyproject.toml migration details (hatchling -> setuptools)
- Test fixture design and golden file format
- Exact CompQuadForm version to extract qfc.c from

</decisions>

<specifics>
## Specific Ideas

- The three-tier fallback chain must match modern tools: GENESIS (TOPMed), GMMAT/SMMAT — these are the gold standard implementations as of 2026
- Davies false convergence is the most dangerous failure mode (ifault=0 but wrong p-value) — the GENESIS+GMMAT hybrid detection strategy addresses this specifically
- The saddlepoint approximation is ~30 lines of math using the full CGF (not moment-matching): K(t) = -1/2 * sum log(1 - 2t*lambda_j), then Lugannani-Rice formula with Newton-Raphson for the saddlepoint equation K'(t_hat) = q
- Key mathematical insight: Liu fails in the tail because no finite set of moments encodes lambda_1 alone; the saddlepoint uses the CGF's singularity at t -> 1/(2*lambda_1) to capture the correct exponential tail rate
- No existing pip-installable pure Python SKAT exists — this would be the first complete implementation with validated accuracy against R
- Cross-language validation literature (pyComBat vs R ComBat) reports ~1e-7 relative squared error; we should achieve better since we use the same qfc.c binary

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 21-pure-python-skat-backend*
*Context gathered: 2026-02-20*
