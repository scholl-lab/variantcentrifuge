# Phase 20: R SKAT Backend - Context

**Gathered:** 2026-02-20
**Status:** Ready for planning

<domain>
## Phase Boundary

Provide SKAT and SKAT-O statistical tests via R's SKAT package through rpy2, with SKATBinary used automatically for binary traits, moment adjustment for small samples, R memory managed explicitly across thousands of genes, and graceful fallback when R is unavailable. The pure Python SKAT backend is Phase 21; ACAT-O combination and diagnostics are Phase 22.

</domain>

<decisions>
## Implementation Decisions

### Fallback behavior
- When R or SKAT is missing: auto-fallback to Python backend with a warning. If Python backend not yet available (Phase 21 not shipped), error with clear message
- When R is installed but SKAT package is missing: same fallback action, but specific message: "R found but SKAT package missing. Install with `install.packages('SKAT')`"
- When R_HOME is misconfigured: show diagnostic error with fix hint (R_HOME value, PATH suggestion), then auto-fallback to Python backend
- Per-gene R crash/exception mid-run: abort entire run immediately — R infrastructure failure is not recoverable
- Statistical NA (too few variants, singular kernel matrix): p_value=NA, log warning, continue — this is a valid statistical outcome, not a crash
- R detection happens eagerly at pipeline startup (when parsing `--skat-backend r`), not lazily at stage execution — fail fast before pipeline work starts

### Output columns & reporting
- `skat_o_rho` column in output — reports the optimal rho value per gene so users can see whether signal is burden-like (rho=1) or variance-like (rho=0)
- `skat_warnings` column capturing R-side warnings per gene (e.g., moment adjustment applied) — useful for auditing
- Aggregate-only timing: log total genes processed and total time at end of R SKAT run. No per-gene timing in normal output
- p_method column design: Claude's Discretion — researcher should investigate SKAT/SAIGE-GENE output conventions and align with existing project column patterns

### Memory & GC tuning
- GC interval configurability: Claude's Discretion — researcher should investigate R/rpy2 memory management best practices and recommend fixed vs configurable, fitting project's existing CLI pattern
- Emit warning before starting R SKAT if gene panel is unusually large (threshold TBD by research, e.g., >2000 genes): "Large gene panel (N genes) — R SKAT memory usage may be significant"
- Monitor R heap size periodically during run. Log WARNING if it exceeds a threshold. Don't abort — just inform
- Progress logging at INFO level every N genes: "SKAT progress: 150/500 genes (30%)" — important for long-running HPC jobs

### R environment detection
- Check R version and SKAT package version at startup. Warn if below recommended minimums (specifics TBD by research) but still proceed
- Comprehensive environment logging when using `--skat-backend r`: R version, SKAT version, rpy2 version, R_HOME path — all at INFO level for reproducibility and HPC debugging

### Claude's Discretion
- p_method column format (per-test vs single backend column) — research best practices from SKAT/SAIGE-GENE output formats
- GC interval: fixed vs configurable — research R/rpy2 memory patterns and recommend
- Minimum R and SKAT version thresholds
- R heap monitoring implementation (what rpy2 exposes, threshold values)
- Progress logging interval (every 10 genes? 50? percentage-based?)

</decisions>

<specifics>
## Specific Ideas

- Fallback philosophy: always try to continue running (fallback to Python, skip statistical NAs) but abort on actual R infrastructure failures — "soft on statistics, hard on crashes"
- HPC reproducibility is a priority — comprehensive R environment logging helps debug issues across cluster nodes with different R installations
- Memory warnings serve as guidance, not gatekeeping — warn but don't block

</specifics>

<deferred>
## Deferred Ideas

None — discussion stayed within phase scope

</deferred>

---

*Phase: 20-r-skat-backend*
*Context gathered: 2026-02-20*
