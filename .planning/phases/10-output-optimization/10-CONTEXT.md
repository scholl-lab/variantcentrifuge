# Phase 10: Output Optimization - Context

**Gathered:** 2026-02-14
**Status:** Ready for planning

<domain>
## Phase Boundary

2-5x faster Excel generation and elimination of redundant GT column parsing across the pipeline. Output must remain functionally equivalent — same data, working links, filters, freeze panes. No new output features or formats.

</domain>

<decisions>
## Implementation Decisions

### Excel output fidelity
- Functional equivalence is sufficient — same data, working hyperlinks, freeze panes, auto-filters
- Minor visual differences from xlsxwriter vs openpyxl (font rendering, column widths) are acceptable
- Excel files are consumed by humans only (clinicians/researchers in Excel or LibreOffice), no programmatic downstream consumers
- Sheet order must be preserved: Results, Metadata, Statistics, Gene Burden
- Memory spike from building all sheets at once is acceptable — optimize for speed over memory

### GT pre-parsed data shape
- Internal cache only — final TSV/Excel output keeps the current GT string format exactly (`Sample(0/1);Sample2(1/1)`)
- Full format fields preserved in cache (DP, AD, etc.), not just genotype call
- Speed priority for storage format — use whatever gives fastest access (dict column, NumPy structured array, etc.)
- Pre-parsing happens at DataFrame load time so all downstream stages benefit

### Hyperlink & formatting scope
- All current features must be preserved: hyperlinks (SpliceAI, Franklin, Varsome, gnomAD, ClinVar), IGV report links, freeze panes, auto-filters
- Link styling doesn't matter as long as clicking opens the right URL — any working link is fine
- No new formatting features (no conditional formatting, no color-coding) — keep output identical to current
- IGV link behavior: just make it work, no specific requirements beyond functional equivalence

### Benchmark strategy
- Use same scales as existing benchmarks: 100, 1K, 10K, 50K variants
- Component + end-to-end benchmarks: measure individual operations (write, hyperlinks, finalize) AND total Excel generation time
- No dedicated GT parsing benchmarks needed — existing stats/inheritance benchmarks prove GT optimization indirectly
- Save benchmark results to .benchmarks/ like Phase 7 and 8 for cross-phase comparison tracking

### Claude's Discretion
- xlsxwriter/openpyxl split strategy (what goes in initial write vs finalization pass)
- GT cache storage format (dict column, structured array, or other)
- Exact benchmark implementation details within the existing framework
- How to batch hyperlink operations for speed
- Column width and styling defaults in xlsxwriter

</decisions>

<specifics>
## Specific Ideas

- Keep output functionally identical — this is a pure speed optimization, not a feature change
- Building all sheets in-memory at once is fine for speed (no incremental append needed)
- GT pre-parsing at load time benefits all downstream stages (stats, inheritance, Excel finalization)
- Follow existing benchmark patterns (100/1K/10K/50K scales, pytest-benchmark framework, saved results)

</specifics>

<deferred>
## Deferred Ideas

- Conditional formatting / color-coding by impact severity — potential future enhancement
- Per-sample GT columns in output — could be useful but changes output format

</deferred>

---

*Phase: 10-output-optimization*
*Context gathered: 2026-02-14*
