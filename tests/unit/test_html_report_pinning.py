"""Unit and regression tests for variant pinning, reviewer notes, and shortlist export (#170)."""

import csv
import http.server
import io
import json
import re
import shutil
import subprocess
import threading
from pathlib import Path

import pytest
from jinja2 import Environment, FileSystemLoader

from variantcentrifuge.generate_html_report import (
    compute_report_id,
    compute_variant_key,
    generate_html_report,
    is_safe_url,
)


def _create_test_env(templates_dir: Path) -> Environment:
    env = Environment(loader=FileSystemLoader(str(templates_dir)))
    env.filters["is_safe_url"] = is_safe_url
    return env


def _run_node_script(js_code: str) -> str:
    """Run a script in node and return stdout."""
    res = subprocess.run(
        ["node"],
        input=js_code,
        capture_output=True,
        text=True,
        check=False,
    )
    if res.returncode != 0:
        raise RuntimeError(
            f"Node script failed (exit {res.returncode}):\n{res.stderr}\nCode:\n{js_code}"
        )
    return res.stdout.strip()


def _dump_dom(chrome_bin: str, url: str) -> subprocess.CompletedProcess:
    """Render ``url`` in headless Chrome and return the final DOM.

    Headless Chrome occasionally hangs in --dump-dom on loaded CI runners; a timeout
    is retried once, while any other failure or assertion still fails immediately.
    """
    cmd = [
        chrome_bin,
        "--headless",
        "--disable-gpu",
        "--no-sandbox",
        "--no-first-run",
        "--virtual-time-budget=6000",
        "--dump-dom",
        url,
    ]
    try:
        return subprocess.run(cmd, capture_output=True, text=True, timeout=45)
    except subprocess.TimeoutExpired:
        return subprocess.run(cmd, capture_output=True, text=True, timeout=45)


def _extract_js_function(html_content: str, fn_name: str) -> str:
    """Extract a JavaScript function by name from rendered HTML report."""
    m = re.search(
        r"(function\s+" + re.escape(fn_name) + r"\s*\([^)]*\)\s*\{[\s\S]*?\n    \})",
        html_content,
    )
    if not m:
        raise ValueError(f"Function {fn_name} not found in rendered HTML")
    return m.group(1)


@pytest.mark.unit
class TestReportIdAndVariantKey:
    """Test deterministic report ID generation and allele-aware variant keys."""

    def test_report_id_deterministic(self) -> None:
        variants = [
            {"CHROM": "1", "POS": 100, "REF": "A", "ALT": "G"},
            {"CHROM": "2", "POS": 200, "REF": "C", "ALT": "T"},
        ]
        id1 = compute_report_id(
            variants, vcf_source="sample.vcf", reference_genome="hg38", filter_expression="None"
        )
        id2 = compute_report_id(
            variants, vcf_source="sample.vcf", reference_genome="hg38", filter_expression="None"
        )
        assert id1 == id2
        assert len(id1) >= 8

        # Explicit config report_id overrides
        id_custom = compute_report_id(
            variants,
            vcf_source="sample.vcf",
            reference_genome="hg38",
            filter_expression="None",
            cfg_report_id="cohort_study_01",
        )
        assert id_custom == "cohort_study_01"

        # Different variants produce different report_id
        variants_diff = [
            {"CHROM": "1", "POS": 100, "REF": "A", "ALT": "C"},
        ]
        id3 = compute_report_id(
            variants_diff,
            vcf_source="sample.vcf",
            reference_genome="hg38",
            filter_expression="None",
        )
        assert id1 != id3

    def test_variant_key_allele_aware(self) -> None:
        # Distinct alternate alleles at same pos have different keys
        v1 = {"CHROM": "chr1", "POS": "1000", "REF": "A", "ALT": "G", "VAR_ID": "rs123"}
        v2 = {"CHROM": "chr1", "POS": "1000", "REF": "A", "ALT": "T", "VAR_ID": "rs123"}
        assert compute_variant_key(v1) != compute_variant_key(v2)

        # Missing or dot VAR_ID does not collide across different loci
        v3 = {"CHROM": "chr1", "POS": "1000", "REF": "A", "ALT": "G", "VAR_ID": "."}
        v4 = {"CHROM": "chr1", "POS": "2000", "REF": "C", "ALT": "T", "VAR_ID": "."}
        assert compute_variant_key(v3) != compute_variant_key(v4)

        # Duplicate annotation rows share the same variant key
        v5 = {"CHROM": "chr1", "POS": "1000", "REF": "A", "ALT": "G", "FEATUREID": "NM_001"}
        v6 = {"CHROM": "chr1", "POS": "1000", "REF": "A", "ALT": "G", "FEATUREID": "NM_001"}
        assert compute_variant_key(v5) == compute_variant_key(v6)

        # Different transcript features for same variant get distinct keys
        v7 = {"CHROM": "chr1", "POS": "1000", "REF": "A", "ALT": "G", "FEATUREID": "NM_002"}
        assert compute_variant_key(v5) != compute_variant_key(v7)

    def test_variant_key_sentinel_normalization_and_hash(self) -> None:
        # Incomplete coordinates with same VAR_ID but different rows produce distinct keys
        v_inc1 = {
            "CHROM": None,
            "POS": ".",
            "REF": "N/A",
            "ALT": "NA",
            "VAR_ID": "rs123",
            "INFO": "a",
        }
        v_inc2 = {
            "CHROM": None,
            "POS": ".",
            "REF": "N/A",
            "ALT": "NA",
            "VAR_ID": "rs123",
            "INFO": "b",
        }
        k1 = compute_variant_key(v_inc1)
        k2 = compute_variant_key(v_inc2)
        assert k1 != k2
        assert k1.startswith("id:rs123:")
        assert k2.startswith("id:rs123:")

        # Incomplete coordinates with missing VAR_ID produce row hash key
        v_inc3 = {"CHROM": None, "POS": ".", "REF": "N/A", "ALT": "NA"}
        assert compute_variant_key(v_inc3).startswith("row:")

        # Sentinel values in feature are normalized to empty feature
        v_norm1 = {"CHROM": "chr1", "POS": "100", "REF": "A", "ALT": "T", "FEATUREID": "."}
        v_norm2 = {"CHROM": "chr1", "POS": "100", "REF": "A", "ALT": "T"}
        assert compute_variant_key(v_norm1) == compute_variant_key(v_norm2) == "chr1:100:A:T"


@pytest.mark.unit
class TestRenderedHTMLPinningElements:
    """Test that generated HTML contains pinning toolbar controls, keys, and embedded records."""

    def test_rendered_html_contains_controls(self, tmp_path: Path) -> None:
        variants = [
            {"CHROM": "17", "POS": "41234567", "REF": "A", "ALT": "G", "GENE": "BRCA1"},
            {"CHROM": "17", "POS": "7577538", "REF": "C", "ALT": "T", "GENE": "TP53"},
        ]
        summary = {"num_variants": 2, "num_genes": 2, "impact_distribution": {}}
        variants_file = tmp_path / "variants.json"
        summary_file = tmp_path / "summary.json"
        variants_file.write_text(json.dumps(variants), encoding="utf-8")
        summary_file.write_text(json.dumps(summary), encoding="utf-8")

        out_dir = tmp_path / "report"
        out_dir.mkdir()
        cfg = {"vcf_file": "/path/to/test.vcf", "reference": "hg38"}

        generate_html_report(str(variants_file), str(summary_file), str(out_dir), cfg)

        html_content = (out_dir / "index.html").read_text(encoding="utf-8")

        # Toolbar controls
        assert 'id="pinned-only-btn"' in html_content
        assert 'id="pinned-count-badge"' in html_content
        assert 'id="export-tsv-btn"' in html_content
        assert 'id="export-json-btn"' in html_content
        assert 'id="clear-review-btn"' in html_content

        # Table rows have data-variant-key
        assert "data-variant-key=" in html_content

        # JavaScript variables
        assert "var reportId =" in html_content
        assert "var rawVariants =" in html_content
        assert "vc_review_" in html_content


@pytest.mark.unit
class TestPinningJSLogicInNode:
    """Test review storage, isolation, escaping, and shortlist generation in node."""

    def test_node_review_state_persistence_and_isolation(self, tmp_path: Path) -> None:
        variants = [{"CHROM": "17", "POS": 41223094, "REF": "T", "ALT": "C", "GENE": "BRCA1"}]
        summary = {"num_variants": 1, "num_genes": 1, "impact_distribution": {}}
        v_file = tmp_path / "v.json"
        s_file = tmp_path / "s.json"
        v_file.write_text(json.dumps(variants), encoding="utf-8")
        s_file.write_text(json.dumps(summary), encoding="utf-8")
        out_dir = tmp_path / "rep"
        out_dir.mkdir()
        generate_html_report(
            str(v_file), str(s_file), str(out_dir), {"vcf_file": "c.vcf", "reference": "hg38"}
        )
        html = (out_dir / "index.html").read_text(encoding="utf-8")

        funcs = "\n\n".join(
            _extract_js_function(html, f)
            for f in [
                "safeStorageGet",
                "safeStorageSet",
                "safeStorageRemove",
                "loadReviewState",
                "saveReviewState",
                "isVariantPinned",
                "getVariantNotes",
                "setVariantPinned",
                "setVariantNotes",
            ]
        )

        script = f"""
        const localStorageStore = {{}};
        const localStorage = {{
            getItem: (k) => localStorageStore[k] || null,
            setItem: (k, v) => {{ localStorageStore[k] = String(v); }},
            removeItem: (k) => {{ delete localStorageStore[k]; }}
        }};
        const window = {{ localStorage }};

        let storageKey = 'vc_review_repA';
        let reviewState = {{ version: 1, variants: {{}} }};
        function updatePinnedUi() {{}}

        {funcs}

        // 1. Initial load
        loadReviewState();
        if (isVariantPinned('17:41223094:T:C')) throw new Error('Initially pinned');

        // 2. Set pinned and notes with production functions
        setVariantPinned('17:41223094:T:C', true);
        setVariantNotes('17:41223094:T:C', 'Pathogenic in BRCA1');

        if (!isVariantPinned('17:41223094:T:C')) {{
            throw new Error('Not pinned after setVariantPinned');
        }}
        if (getVariantNotes('17:41223094:T:C') !== 'Pathogenic in BRCA1') {{
            throw new Error('Notes not set');
        }}

        // 3. Test report isolation: Report B has no review state
        storageKey = 'vc_review_repB';
        loadReviewState();
        if (isVariantPinned('17:41223094:T:C')) {{
            throw new Error('Report B leaked state from Report A');
        }}
        if (getVariantNotes('17:41223094:T:C') !== '') {{
            throw new Error('Report B leaked notes');
        }}

        // 4. Reload Report A from storage
        storageKey = 'vc_review_repA';
        reviewState = {{ version: 1, variants: {{}} }};
        loadReviewState();
        if (!isVariantPinned('17:41223094:T:C')) {{
            throw new Error('Report A did not persist pinned state');
        }}
        if (getVariantNotes('17:41223094:T:C') !== 'Pathogenic in BRCA1') {{
            throw new Error('Report A did not persist notes');
        }}

        // 5. Clear Report A
        safeStorageRemove('vc_review_repA');
        loadReviewState();
        if (isVariantPinned('17:41223094:T:C')) throw new Error('Report A not cleared');

        console.log('OK');
        """
        assert _run_node_script(script) == "OK"

    def test_node_malformed_and_quota_error_storage_resilience(self, tmp_path: Path) -> None:
        variants = [{"CHROM": "1", "POS": 100, "REF": "A", "ALT": "G"}]
        summary = {"num_variants": 1, "num_genes": 1, "impact_distribution": {}}
        v_file = tmp_path / "v.json"
        s_file = tmp_path / "s.json"
        v_file.write_text(json.dumps(variants), encoding="utf-8")
        s_file.write_text(json.dumps(summary), encoding="utf-8")
        out_dir = tmp_path / "rep"
        out_dir.mkdir()
        generate_html_report(str(v_file), str(s_file), str(out_dir), {})
        html = (out_dir / "index.html").read_text(encoding="utf-8")

        funcs = "\n\n".join(
            _extract_js_function(html, f)
            for f in [
                "safeStorageGet",
                "safeStorageSet",
                "safeStorageRemove",
                "loadReviewState",
                "saveReviewState",
                "isVariantPinned",
                "getVariantNotes",
                "setVariantPinned",
                "setVariantNotes",
            ]
        )

        script = f"""
        let mockStorage = {{}};
        let getItemOverride = null;
        let setItemOverride = null;
        let removeItemOverride = null;

        const localStorage = {{
            getItem: (k) => getItemOverride ? getItemOverride(k) : (mockStorage[k] || null),
            setItem: (k, v) => {{
                if (setItemOverride) return setItemOverride(k, v);
                mockStorage[k] = String(v);
            }},
            removeItem: (k) => {{
                if (removeItemOverride) return removeItemOverride(k);
                delete mockStorage[k];
            }}
        }};
        const window = {{ localStorage }};

        let storageKey = 'vc_review_malformed';
        let reviewState = {{ version: 1, variants: {{}} }};
        function updatePinnedUi() {{}}

        {funcs}

        // 1. Corrupt JSON string recovery
        mockStorage[storageKey] = '{{corrupt_json:';
        loadReviewState();
        if (typeof reviewState.variants !== 'object' || Array.isArray(reviewState.variants)) {{
            throw new Error('Did not recover from corrupt JSON');
        }}

        // 2. Malformed variants array: {{\"version\": 1, \"variants\": []}}
        mockStorage[storageKey] = JSON.stringify({{ version: 1, variants: [] }});
        loadReviewState();
        if (typeof reviewState.variants !== 'object' || Array.isArray(reviewState.variants)) {{
            throw new Error('Did not handle variants as array');
        }}

        // 3. Malformed entries inside variants object
        mockStorage[storageKey] = JSON.stringify({{
            version: 1,
            variants: {{
                bad_array: [1, 2],
                bad_null: null,
                bad_str: 'just text',
                good_entry: {{ pinned: 1, notes: 999 }}
            }}
        }});
        loadReviewState();
        var bad1 = reviewState.variants['bad_array'];
        var bad2 = reviewState.variants['bad_null'];
        var bad3 = reviewState.variants['bad_str'];
        if (bad1 || bad2 || bad3) {{
            throw new Error('Malformed variant entries were not rejected');
        }}
        var good = reviewState.variants['good_entry'];
        if (!good || good.pinned !== true || good.notes !== '') {{
            throw new Error('Sanitization failed on good_entry');
        }}

        // 4. QuotaExceededError on save: safeStorageSet returns false, note preserved in session
        setItemOverride = function(k, v) {{
            const err = new Error('Quota exceeded');
            err.name = 'QuotaExceededError';
            throw err;
        }};

        const persistResult = setVariantNotes('good_entry', 'offline note');
        if (persistResult !== false) {{
            throw new Error('Expected persistResult to be false on QuotaExceededError');
        }}
        if (getVariantNotes('good_entry') !== 'offline note') {{
            throw new Error('In-memory notes were lost on QuotaExceededError');
        }}

        // 5. SecurityError / Storage disabled: safe storage helpers return null/false
        const secErr = () => {{
            const err = new Error('Access denied');
            err.name = 'SecurityError';
            throw err;
        }};
        getItemOverride = secErr;
        setItemOverride = secErr;
        removeItemOverride = secErr;

        if (safeStorageGet('key') !== null) {{
            throw new Error('safeStorageGet threw or did not return null');
        }}
        if (safeStorageSet('key', 'val') !== false) {{
            throw new Error('safeStorageSet threw or did not return false');
        }}
        if (safeStorageRemove('key') !== false) {{
            throw new Error('safeStorageRemove threw or did not return false');
        }}

        console.log('OK');
        """
        assert _run_node_script(script) == "OK"

    def test_export_shortlist_tsv_and_json_serialization(self, tmp_path: Path) -> None:
        """Verify production TSV and JSON export logic with special characters and Unicode."""
        variants = [
            {
                "CHROM": "chr17",
                "POS": 41234567,
                "REF": "A",
                "ALT": "G",
                "GENE": "BRCA1",
                "nested_meta": {"score": 0.99, "call": "PASS"},
                "igv_links": [{"report_path": "igv.html"}],
            },
            {
                "CHROM": "chr7",
                "POS": 140453136,
                "REF": "A",
                "ALT": "T",
                "GENE": "BRAF",
                "igv_links": [],
            },
            {
                "CHROM": "chr13",
                "POS": 32900000,
                "REF": "C",
                "ALT": "T",
                "GENE": "BRCA2",
                "igv_links": [],
            },
        ]
        summary = {"num_variants": 3, "num_genes": 3, "impact_distribution": {}}
        v_file = tmp_path / "v.json"
        s_file = tmp_path / "s.json"
        v_file.write_text(json.dumps(variants), encoding="utf-8")
        s_file.write_text(json.dumps(summary), encoding="utf-8")
        out_dir = tmp_path / "rep"
        out_dir.mkdir()
        generate_html_report(
            str(v_file),
            str(s_file),
            str(out_dir),
            {"vcf_file": "cohort.vcf", "reference": "hg38"},
        )
        html = (out_dir / "index.html").read_text(encoding="utf-8")

        funcs = "\n\n".join(
            _extract_js_function(html, f)
            for f in [
                "isVariantPinned",
                "escapeTsvField",
                "getShortlistedCandidates",
                "generateShortlistTsv",
                "generateShortlistJson",
            ]
        )

        m_var = re.search(r"var rawVariants = (\[[\s\S]*?\]);", html)
        m_rep = re.search(r'var reportId = ("[^"]*");', html)
        m_cols = re.search(r"var allColumnOriginalNames = (\[[^\n]*\]);", html)
        assert m_var and m_rep and m_cols

        note_text = (
            'Validated pathogenic: "c.123A>G"\t(ClinVar)\nNote: \u03b1/\u03b2 kinase '
            "& <script>alert(1)</script>"
        )

        script = f"""
        const rawVariants = {m_var.group(1)};
        const reportId = {m_rep.group(1)};
        const allColumnOriginalNames = {m_cols.group(1)};
        const reviewState = {{
            version: 1,
            variants: {{
                'chr17:41234567:A:G': {{
                    pinned: true,
                    notes: {json.dumps(note_text)},
                    updated_at: '2026-10-19T20:00:00.000Z'
                }},
                'chr7:140453136:A:T': {{
                    pinned: true,
                    notes: 'V600E candidate',
                    updated_at: '2026-10-19T20:05:00.000Z'
                }}
            }}
        }};

        {funcs}

        const tsv = generateShortlistTsv();
        const jsonStr = generateShortlistJson();

        console.log(JSON.stringify({{ tsv, json: jsonStr }}));
        """

        output_json = _run_node_script(script)
        data = json.loads(output_json)
        tsv_content = data["tsv"]
        json_content = data["json"]

        # Assert no internal keys or [object Object]
        assert "_variant_key" not in tsv_content
        assert "igv_links" not in tsv_content
        assert "[object Object]" not in tsv_content

        # 1. Verify TSV parses correctly using Python standard csv reader
        reader = list(csv.reader(io.StringIO(tsv_content), delimiter="\t"))
        assert len(reader) == 3  # Header + 2 candidates (chr13 unpinned is excluded)
        header = reader[0]
        assert "CHROM" in header
        assert "nested_meta" in header
        assert "review_notes" in header
        assert "review_timestamp" in header
        assert "report_id" in header

        row1 = dict(zip(header, reader[1], strict=False))
        assert row1["CHROM"] == "chr17"
        assert row1["GENE"] == "BRCA1"
        assert json.loads(row1["nested_meta"]) == {"score": 0.99, "call": "PASS"}
        assert "\t(ClinVar)" in row1["review_notes"]
        assert '"c.123A>G"' in row1["review_notes"]
        assert "\nNote:" in row1["review_notes"]
        assert "\u03b1/\u03b2 kinase" in row1["review_notes"]
        assert "<script>alert(1)</script>" in row1["review_notes"]

        # 2. Verify JSON format
        parsed_json = json.loads(json_content)
        assert parsed_json["candidate_count"] == 2
        candidates = parsed_json["candidates"]
        assert len(candidates) == 2
        assert "_variant_key" not in candidates[0]
        assert "igv_links" not in candidates[0]
        assert candidates[0]["GENE"] == "BRCA1"
        assert candidates[0]["nested_meta"] == {"score": 0.99, "call": "PASS"}
        assert "\u03b1/\u03b2 kinase" in candidates[0]["review_notes"]
        assert "<script>" in candidates[0]["review_notes"]

    def test_notes_containing_markup_remain_inert_in_rendered_html(self, tmp_path: Path) -> None:
        """Verify notes with html/script markup are escaped and remain inert text in child rows."""
        variants = [{"CHROM": "1", "POS": "100", "REF": "A", "ALT": "G", "GENE": "TEST"}]
        summary = {"num_variants": 1, "num_genes": 1, "impact_distribution": {}}
        v_file = tmp_path / "v.json"
        s_file = tmp_path / "s.json"
        v_file.write_text(json.dumps(variants), encoding="utf-8")
        s_file.write_text(json.dumps(summary), encoding="utf-8")
        out_dir = tmp_path / "rep"
        out_dir.mkdir()
        generate_html_report(str(v_file), str(s_file), str(out_dir), {})
        html = (out_dir / "index.html").read_text(encoding="utf-8")

        script = f"""
        function escapeHtml(str) {{
            if (str === null || str === undefined) return '';
            return String(str)
                .replace(/&/g, '&amp;')
                .replace(/</g, '&lt;')
                .replace(/>/g, '&gt;')
                .replace(/"/g, '&quot;')
                .replace(/'/g, '&#39;');
        }}
        function isSafeUrl(u) {{ return true; }}
        function extractCellContent(d) {{ return {{ text: String(d || ''), links: [] }}; }}

        const reviewState = {{
            version: 1,
            variants: {{
                "1:100:A:G": {{
                    pinned: true,
                    notes: '<script>alert("xss")</script><img src=x onerror=alert(2)>',
                    updated_at: '2026-10-19T20:00:00Z'
                }}
            }}
        }};
        function isVariantPinned(k) {{
            return Boolean(reviewState.variants[k] && reviewState.variants[k].pinned);
        }}
        function getVariantNotes(k) {{
            return (reviewState.variants[k] && reviewState.variants[k].notes) || '';
        }}

        // Extract formatChildRow from rendered html
        {html[html.find("function formatChildRow(") : html.find("// Initialize DataTables v2")]}

        const rowData = [null, 'TEST', '1', '100', 'A', 'G'];
        const columns = [
            {{ original_name: 'GENE' }},
            {{ original_name: 'CHROM' }},
            {{ original_name: 'POS' }},
            {{ original_name: 'REF' }},
            {{ original_name: 'ALT' }}
        ];
        const res = formatChildRow(rowData, columns, "1:100:A:G");

        if (res.includes('<script>') || res.includes('<img')) {{
            console.error('Markup in note was not neutralized:', res);
            process.exit(1);
        }}
        if (!res.includes('&lt;script&gt;alert(&quot;xss&quot;)&lt;/script&gt;')) {{
            console.error('Note content missing or improperly escaped:', res);
            process.exit(1);
        }}
        console.log('OK');
        """
        assert _run_node_script(script) == "OK"


@pytest.mark.unit
class TestHeadlessChromeInteractivePinning:
    """Browser-interactive test for pinning, notes, filtering, and export in headless Chrome."""

    def test_interactive_pinning_and_shortlist_export(self, tmp_path: Path) -> None:
        chrome_bin = (
            shutil.which("google-chrome")
            or shutil.which("chromium")
            or shutil.which("chromium-browser")
        )
        if not chrome_bin:
            pytest.skip("Chrome/Chromium binary not available")

        # 16 variants across 2 pages with duplicate annotation rows across pages
        # and sentinel coordinates
        variants_data = [
            {
                "CHROM": "chr17",
                "POS": 41223094,
                "REF": "T",
                "ALT": "C",
                "GENE": "BRCA1",
                "IMPACT": "HIGH",
                "VAR_ID": "rs123",
            },
            {
                "CHROM": "chr13",
                "POS": 32900000,
                "REF": "C",
                "ALT": "T",
                "GENE": "BRCA2",
                "FEATUREID": "NM_000059",
                "IMPACT": "MODERATE",
            },
            {
                "CHROM": "chr13",
                "POS": 32900000,
                "REF": "C",
                "ALT": "T",
                "GENE": "BRCA2",
                "FEATUREID": "NM_000059",
                "IMPACT": "LOW",
            },
            {
                "CHROM": "chr1",
                "POS": 1000,
                "REF": "A",
                "ALT": "T",
                "GENE": "GENE1",
                "IMPACT": "MODERATE",
            },
            {
                "CHROM": "chr1",
                "POS": 1000,
                "REF": "A",
                "ALT": "C",
                "GENE": "GENE1",
                "IMPACT": "MODERATE",
            },
            {
                "CHROM": None,
                "POS": ".",
                "REF": "N/A",
                "ALT": "NA",
                "GENE": "GENE2",
                "VAR_ID": "rs999",
                "EXTRA": "rowA",
            },
            {
                "CHROM": None,
                "POS": ".",
                "REF": "N/A",
                "ALT": "NA",
                "GENE": "GENE3",
                "VAR_ID": "rs999",
                "EXTRA": "rowB",
            },
            {
                "CHROM": "chr2",
                "POS": 2000,
                "REF": "A",
                "ALT": "G",
                "GENE": "GENE4",
                "IMPACT": "LOW",
            },
            {
                "CHROM": "chr3",
                "POS": 3000,
                "REF": "A",
                "ALT": "G",
                "GENE": "GENE5",
                "IMPACT": "LOW",
            },
            {
                "CHROM": "chr4",
                "POS": 4000,
                "REF": "A",
                "ALT": "G",
                "GENE": "GENE6",
                "IMPACT": "LOW",
            },
            {
                "CHROM": "chr13",
                "POS": 32900000,
                "REF": "C",
                "ALT": "T",
                "GENE": "BRCA2",
                "FEATUREID": "NM_000059",
                "IMPACT": "MODIFIER",
            },
            {
                "CHROM": "chr17",
                "POS": 7577538,
                "REF": "C",
                "ALT": "T",
                "GENE": "TP53",
                "IMPACT": "HIGH",
            },
            {
                "CHROM": "chr5",
                "POS": 5000,
                "REF": "A",
                "ALT": "G",
                "GENE": "GENE7",
                "IMPACT": "LOW",
            },
            {
                "CHROM": "chr6",
                "POS": 6000,
                "REF": "A",
                "ALT": "G",
                "GENE": "GENE8",
                "IMPACT": "LOW",
            },
            {
                "CHROM": "chr8",
                "POS": 8000,
                "REF": "A",
                "ALT": "G",
                "GENE": "GENE9",
                "IMPACT": "LOW",
            },
            {
                "CHROM": "chr9",
                "POS": 9000,
                "REF": "A",
                "ALT": "G",
                "GENE": "GENE10",
                "IMPACT": "LOW",
            },
        ]
        summary = {
            "num_variants": 16,
            "num_genes": 11,
            "num_samples": 1,
            "impact_distribution": {"HIGH": 2, "MODERATE": 3, "LOW": 10, "MODIFIER": 1},
        }
        cfg = {"vcf_file": "cohort.vcf", "reference": "hg38"}
        v_file = tmp_path / "variants.json"
        s_file = tmp_path / "summary.json"
        v_file.write_text(json.dumps(variants_data), encoding="utf-8")
        s_file.write_text(json.dumps(summary), encoding="utf-8")

        generate_html_report(str(v_file), str(s_file), str(tmp_path), cfg)
        report_html = (tmp_path / "index.html").read_text(encoding="utf-8")

        probe_script = """
<script>
window.addEventListener("load", function() {
    setTimeout(function() {
        try {
            var dt = $("#variants_table").DataTable();
            // Set page length to 10 so we test multi-page pagination
            dt.page.len(10).draw();
            var totalRows = dt.rows().count();
            var countBadge = document.getElementById("pinned-count-badge");
            var initBadge = countBadge ? countBadge.textContent.trim() : "";

            // 1. Keyboard accessibility: Toggle row details using Enter and Space keys
            var trs = $("#variants_table tbody tr[data-variant-key]");
            var row0Tr = $(trs[0]);
            var row0ControlTd = row0Tr.find("td.dt-control, td.details-control");
            var row0ToggleBtn = row0Tr.find(".details-toggle");

            // (a) Enter key on td.dt-control expands child row
            row0ControlTd[0].dispatchEvent(
                new KeyboardEvent("keydown", { key: "Enter", keyCode: 13, bubbles: true })
            );
            var expandedAfterEnter = row0Tr.hasClass("shown") && (
                row0ControlTd.attr("aria-expanded") === "true" ||
                row0ToggleBtn.attr("aria-expanded") === "true"
            );

            // (b) Space key on td.dt-control collapses child row
            row0ControlTd[0].dispatchEvent(
                new KeyboardEvent("keydown", { key: " ", keyCode: 32, bubbles: true })
            );
            var collapsedAfterSpace = !row0Tr.hasClass("shown") && (
                row0ControlTd.attr("aria-expanded") === "false" ||
                row0ToggleBtn.attr("aria-expanded") === "false"
            );

            // (c) Enter key on .details-toggle button expands child row again
            row0ToggleBtn[0].dispatchEvent(
                new KeyboardEvent("keydown", { key: "Enter", keyCode: 13, bubbles: true })
            );
            var expandedAfterBtnEnter = row0Tr.hasClass("shown") && (
                row0ControlTd.attr("aria-expanded") === "true" ||
                row0ToggleBtn.attr("aria-expanded") === "true"
            );

            // Initial pin and notes on page 1 (row 0: BRCA1)
            var pinBtn0 = row0Tr.next().find(".btn-pin-variant, .pin-toggle-btn");
            pinBtn0.trigger("click");

            var noteInput0 = row0Tr.next().find(".review-notes-input");
            var inertMarkupText = (
                '<script>alert("xss")<\\/script><b>Pathogenic</b>\\t"c.123A>G"\\n' +
                "Note: \\u03b1/\\u03b2 kinase"
            );
            noteInput0.val(inertMarkupText);
            noteInput0[0].dispatchEvent(new Event("input", { bubbles: true }));

            var sec0 = noteInput0.closest(".review-section");
            var status0 = sec0.find(".review-status").text();

            // 2. Off-page duplicate child row synchronization across pages:
            // First, switch to page 2 (indices 10..15) and expand duplicate row 10 (BRCA2)
            dt.page(1).draw("page");
            var page2Trs = $("#variants_table tbody tr[data-variant-key]");
            var p2Brca2Tr = null;
            page2Trs.each(function() {
                if ($(this).attr("data-variant-key") === "chr13:32900000:C:T:NM_000059") {
                    p2Brca2Tr = $(this);
                }
            });
            if (p2Brca2Tr) {
                p2Brca2Tr.find("td.dt-control, td.details-control").trigger("click");
            }

            // Return to page 1: child row for page 2 BRCA2 is now cached in DT memory!
            dt.page(0).draw("page");
            var p1Trs = $("#variants_table tbody tr[data-variant-key]");
            var brca2Tr1 = $(p1Trs[1]);
            brca2Tr1.find("td.dt-control, td.details-control").trigger("click");
            var brca2PinBtn1 = brca2Tr1.next().find(".btn-pin-variant, .pin-toggle-btn");
            brca2PinBtn1.trigger("click");

            var brca2NoteInput1 = brca2Tr1.next().find(".review-notes-input");
            brca2NoteInput1.val("Duplicate across page 1 and page 2");
            brca2NoteInput1[0].dispatchEvent(new Event("input", { bubbles: true }));

            // Open duplicate row 2 on page 1: verify pin button and note value synced in DOM
            var brca2Tr2 = $(p1Trs[2]);
            brca2Tr2.find("td.dt-control, td.details-control").trigger("click");
            var brca2PinBtn2 = brca2Tr2.next().find(".btn-pin-variant, .pin-toggle-btn");
            var brca2NoteInput2 = brca2Tr2.next().find(".review-notes-input");
            var isDupBtnPinned = brca2PinBtn2.hasClass("pinned");
            var dupNoteVal = brca2NoteInput2.val();

            // Verify detached/cached DataTables child row for the off-page duplicate row (index 10)
            var p2Row = dt.row(10);
            var p2Child = p2Row.child();
            var cachedP2BtnPinned = false;
            var cachedP2NoteVal = "";
            if (p2Child && p2Child.length) {
                var cachedBtn = $(p2Child[0]).find(".btn-pin-variant, .pin-toggle-btn");
                cachedP2BtnPinned = (
                    cachedBtn.hasClass("pinned") && cachedBtn.attr("aria-pressed") === "true"
                );
                cachedP2NoteVal = $(p2Child[0]).find(".review-notes-input").val();
            }

            // Now navigate to page 2 and verify live DOM on page 2:
            dt.page(1).draw("page");
            var liveP2Tr = null;
            $("#variants_table tbody tr[data-variant-key]").each(function() {
                if ($(this).attr("data-variant-key") === "chr13:32900000:C:T:NM_000059") {
                    liveP2Tr = $(this);
                }
            });
            var liveP2BtnPinned = false;
            var liveP2NoteVal = "";
            if (liveP2Tr) {
                var liveBtn = liveP2Tr.next().find(".btn-pin-variant, .pin-toggle-btn");
                liveP2BtnPinned = (
                    liveBtn.hasClass("pinned") && liveBtn.attr("aria-pressed") === "true"
                );
                liveP2NoteVal = liveP2Tr.next().find(".review-notes-input").val();
            }

            // Also pin TP53 on page 2
            var tp53Tr = null;
            $("#variants_table tbody tr[data-variant-key]").each(function() {
                if ($(this).text().indexOf("TP53") !== -1) tp53Tr = $(this);
            });
            if (tp53Tr) {
                tp53Tr.find("td.dt-control, td.details-control").trigger("click");
                var tp53PinBtn = tp53Tr.next().find(".btn-pin-variant, .pin-toggle-btn");
                tp53PinBtn.trigger("click");
                var tp53Note = tp53Tr.next().find(".review-notes-input");
                tp53Note.val("TP53 variant on page 2");
                tp53Note[0].dispatchEvent(new Event("input", { bubbles: true }));
            }
            var badgeAfterPage2 = countBadge ? countBadge.textContent.trim() : "";

            // Return to page 1
            dt.page(0).draw("page");

            // 3. Repeated VAR_ID with missing coords: row 5 and row 6
            var repTrs = $("#variants_table tbody tr[data-variant-key]");
            var repIdTr5 = $(repTrs[5]);
            var repIdTr6 = $(repTrs[6]);
            var key5 = repIdTr5.attr("data-variant-key");
            var key6 = repIdTr6.attr("data-variant-key");
            var keysDistinct = (
                (key5 !== key6) &&
                (key5.indexOf("id:rs999:") === 0) &&
                (key6.indexOf("id:rs999:") === 0)
            );

            // 4. Distinct alternate alleles: row 3 and row 4
            var altTr3 = $(repTrs[3]);
            var altTr4 = $(repTrs[4]);
            var key3 = altTr3.attr("data-variant-key");
            var key4 = altTr4.attr("data-variant-key");
            var altKeysDistinct = (key3 !== key4);

            // 5. Column visibility change and sorting survival
            dt.order([[1, "desc"]]).draw();
            dt.column(2).visible(false);

            // 6. Search filter 'BRCA' (TP53 on page 2 filtered out)
            dt.search("BRCA").draw();
            var searchAppliedCount = dt.rows({ search: "applied" }).count();

            // 7. TSV and JSON exports contain both visible and filtered-out/off-page candidates
            var tsv = window.generateShortlistTsv();
            var jsonStr = window.generateShortlistJson();

            // 8. Storage resilience: QuotaExceededError
            var origSetItem = window.localStorage.setItem;
            window.localStorage.setItem = function() {
                var err = new Error("Quota exceeded");
                err.name = "QuotaExceededError";
                throw err;
            };
            var testInput = $(".review-notes-input").first();
            if (testInput.length > 0) {
                testInput.val("Offline session note");
                testInput[0].dispatchEvent(new Event("input", { bubbles: true }));
            }
            var offlineStatus = testInput.closest(".review-section").find(".review-status").text();
            window.localStorage.setItem = origSetItem;

            // 9. Report isolation: external report review key does not affect this report
            var unrelatedData = JSON.stringify({
                version: 1,
                variants: { "other:1:A:G": { pinned: true, notes: "isolated" } }
            });
            window.safeStorageSet("vc_review_unrelated_report_xyz", unrelatedData);
            window.loadReviewState();
            var badgeAfterIso = countBadge ? countBadge.textContent.trim() : "";

            // Reset search filter so all rows are visible
            dt.search("").draw();

            // 10. Clear Review with populated child rows on both page 1 and page 2:
            // Switch to page 2 (BRCA2 and TP53 rows are open)
            dt.page(1).draw("page");
            window.confirm = function() { return true; };
            var clearBtn = document.getElementById("clear-review-btn");
            if (clearBtn) clearBtn.click();
            var badgeAfterClear = countBadge ? countBadge.textContent.trim() : "";

            // Check page 2 child rows are cleared of pin and notes:
            var pinnedSel = (
                "#variants_table tbody .btn-pin-variant.pinned, " +
                "#variants_table tbody .pin-toggle-btn.pinned"
            );
            var page2BtnsStillPinned = $(pinnedSel).length;
            var page2NotesCleared = true;
            $("#variants_table tbody .review-notes-input").each(function() {
                if ($(this).val() !== "") page2NotesCleared = false;
            });

            // Switch to page 1: check page 1 child rows are also cleared:
            dt.page(0).draw("page");
            var page1BtnsStillPinned = $(pinnedSel).length;
            var page1NotesCleared = true;
            $("#variants_table tbody .review-notes-input").each(function() {
                if ($(this).val() !== "") page1NotesCleared = false;
            });

            var result = {
                totalRows: totalRows,
                initBadge: initBadge,
                expandedAfterEnter: expandedAfterEnter,
                collapsedAfterSpace: collapsedAfterSpace,
                expandedAfterBtnEnter: expandedAfterBtnEnter,
                status0: status0,
                isDupBtnPinned: isDupBtnPinned,
                dupNoteVal: dupNoteVal,
                cachedP2BtnPinned: cachedP2BtnPinned,
                cachedP2NoteVal: cachedP2NoteVal,
                liveP2BtnPinned: liveP2BtnPinned,
                liveP2NoteVal: liveP2NoteVal,
                keysDistinct: keysDistinct,
                altKeysDistinct: altKeysDistinct,
                badgeAfterPage2: badgeAfterPage2,
                searchAppliedCount: searchAppliedCount,
                hasBrcaInTsv: tsv.indexOf("BRCA1") !== -1,
                hasTp53InTsv: tsv.indexOf("TP53") !== -1,
                hasBrcaInJson: jsonStr.indexOf("BRCA1") !== -1,
                hasTp53InJson: jsonStr.indexOf("TP53") !== -1,
                hasObjectObjectInTsv: tsv.indexOf("[object Object]") !== -1,
                hasVariantKeyInTsv: tsv.indexOf("_variant_key") !== -1,
                offlineStatus: offlineStatus,
                badgeAfterIso: badgeAfterIso,
                badgeAfterClear: badgeAfterClear,
                page2BtnsStillPinned: page2BtnsStillPinned,
                page2NotesCleared: page2NotesCleared,
                page1BtnsStillPinned: page1BtnsStillPinned,
                page1NotesCleared: page1NotesCleared
            };

            var marker = document.createElement("div");
            marker.id = "browser-probe-result";
            marker.textContent = JSON.stringify(result);
            document.body.appendChild(marker);
        } catch (e) {
            var errMarker = document.createElement("div");
            errMarker.id = "browser-probe-result";
            errMarker.textContent = JSON.stringify({ error: e.message, stack: e.stack });
            document.body.appendChild(errMarker);
        }
    }, 350);
});
</script>
"""
        probed_html = report_html.replace("</body>", probe_script + "</body>")

        class Handler(http.server.BaseHTTPRequestHandler):
            def do_GET(self):
                self.send_response(200)
                self.send_header("Content-Type", "text/html; charset=utf-8")
                self.end_headers()
                self.wfile.write(probed_html.encode("utf-8"))

            def log_message(self, *args):
                pass

        server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Handler)
        port = server.server_port
        t = threading.Thread(target=server.serve_forever)
        t.daemon = True
        t.start()

        try:
            res = _dump_dom(chrome_bin, f"http://127.0.0.1:{port}/")
            assert res.returncode == 0, f"Chrome failed: {res.stderr}"
            m = re.search(r'<div id="browser-probe-result">([\s\S]*?)</div>', res.stdout)
            assert m is not None, f"Marker not found: {res.stdout[:500]}"
            import html as html_lib

            data = json.loads(html_lib.unescape(m.group(1)))
            assert "error" not in data, f"Browser error: {data.get('error')}\n{data.get('stack')}"

            assert data["totalRows"] == 16
            assert data["initBadge"] == "0"
            # Keyboard row toggle assertions
            assert data["expandedAfterEnter"] is True
            assert data["collapsedAfterSpace"] is True
            assert data["expandedAfterBtnEnter"] is True
            assert data["status0"] == "Saved"

            # Duplicate row sync assertions (page 1 and cached page 2)
            assert data["isDupBtnPinned"] is True
            assert data["dupNoteVal"] == "Duplicate across page 1 and page 2"
            assert data["cachedP2BtnPinned"] is True
            assert data["cachedP2NoteVal"] == "Duplicate across page 1 and page 2"
            assert data["liveP2BtnPinned"] is True
            assert data["liveP2NoteVal"] == "Duplicate across page 1 and page 2"

            assert data["keysDistinct"] is True
            assert data["altKeysDistinct"] is True
            assert data["badgeAfterPage2"] == "3"
            assert data["searchAppliedCount"] < 16

            # TSV and JSON exports contain BOTH visible (BRCA1) and off-page/filtered-out (TP53)
            assert data["hasBrcaInTsv"] is True
            assert data["hasTp53InTsv"] is True
            assert data["hasBrcaInJson"] is True
            assert data["hasTp53InJson"] is True
            assert data["hasObjectObjectInTsv"] is False
            assert data["hasVariantKeyInTsv"] is False

            # Offline / quota exceeded status feedback
            assert data["offlineStatus"] == "Saved in session only"
            # Report isolation preserved
            assert data["badgeAfterIso"] == "3"

            # Clear review resets both active and off-page child rows
            assert data["badgeAfterClear"] == "0"
            assert data["page2BtnsStillPinned"] == 0
            assert data["page2NotesCleared"] is True
            assert data["page1BtnsStillPinned"] == 0
            assert data["page1NotesCleared"] is True
        finally:
            server.shutdown()
            server.server_close()

    def test_headless_chrome_reload_and_multireport_isolation(self, tmp_path: Path) -> None:
        chrome_bin = (
            shutil.which("google-chrome")
            or shutil.which("chromium")
            or shutil.which("chromium-browser")
        )
        if not chrome_bin:
            pytest.skip("Chrome/Chromium binary not available")

        # Generate Report 1
        rep1_variants = [
            {
                "CHROM": "17",
                "POS": 41223094,
                "REF": "T",
                "ALT": "C",
                "GENE": "BRCA1",
                "IMPACT": "HIGH",
            },
            {
                "CHROM": "17",
                "POS": 7577538,
                "REF": "C",
                "ALT": "T",
                "GENE": "TP53",
                "IMPACT": "HIGH",
            },
        ]
        rep1_dir = tmp_path / "rep1"
        rep1_dir.mkdir()
        (rep1_dir / "v.json").write_text(json.dumps(rep1_variants), encoding="utf-8")
        (rep1_dir / "s.json").write_text(
            json.dumps({"num_variants": 2, "num_genes": 2, "impact_distribution": {}}),
            encoding="utf-8",
        )
        generate_html_report(
            str(rep1_dir / "v.json"),
            str(rep1_dir / "s.json"),
            str(rep1_dir),
            {"vcf_file": "c1.vcf", "reference": "hg38", "report_id": "report_alpha"},
        )
        rep1_html = (rep1_dir / "index.html").read_text(encoding="utf-8")

        # Generate Report 2
        rep2_variants = [
            {
                "CHROM": "7",
                "POS": 55241707,
                "REF": "G",
                "ALT": "A",
                "GENE": "EGFR",
                "IMPACT": "HIGH",
            },
            {
                "CHROM": "12",
                "POS": 25398284,
                "REF": "C",
                "ALT": "T",
                "GENE": "KRAS",
                "IMPACT": "HIGH",
            },
        ]
        rep2_dir = tmp_path / "rep2"
        rep2_dir.mkdir()
        (rep2_dir / "v.json").write_text(json.dumps(rep2_variants), encoding="utf-8")
        (rep2_dir / "s.json").write_text(
            json.dumps({"num_variants": 2, "num_genes": 2, "impact_distribution": {}}),
            encoding="utf-8",
        )
        generate_html_report(
            str(rep2_dir / "v.json"),
            str(rep2_dir / "s.json"),
            str(rep2_dir),
            {"vcf_file": "c2.vcf", "reference": "hg38", "report_id": "report_beta"},
        )
        rep2_html = (rep2_dir / "index.html").read_text(encoding="utf-8")

        rep1_probe = """
<script>
window.addEventListener("load", function() {
    setTimeout(function() {
        try {
            var step = sessionStorage.getItem("reload_test_step") || "1";
            if (step === "1") {
                // Step 1: Pin BRCA1 and write note on Report 1
                var tr0 = $("#variants_table tbody tr[data-variant-key]").first();
                tr0.find("td.dt-control, td.details-control").trigger("click");
                var pinBtn = tr0.next().find(".btn-pin-variant, .pin-toggle-btn");
                pinBtn.trigger("click");
                var noteInp = tr0.next().find(".review-notes-input");
                noteInp.val("Persistent review note for BRCA1 across reload");
                noteInp[0].dispatchEvent(new Event("input", { bubbles: true }));

                var countBadge = document.getElementById("pinned-count-badge");
                var badgeStep1 = countBadge ? countBadge.textContent.trim() : "";
                sessionStorage.setItem("reload_test_step", "2");
                sessionStorage.setItem("rep1_badge_step1", badgeStep1);

                // Perform real window reload in headless Chrome
                window.location.reload();
                return;
            } else if (step === "2") {
                // Step 2: Running after window.location.reload()
                var countBadge = document.getElementById("pinned-count-badge");
                var badgeAfterReload = countBadge ? countBadge.textContent.trim() : "";
                var isPinnedAfterReload = window.isVariantPinned("17:41223094:T:C");
                var noteAfterReload = window.getVariantNotes("17:41223094:T:C");

                sessionStorage.setItem("reload_test_step", "3");
                sessionStorage.setItem("rep1_badge_after_reload", badgeAfterReload);
                sessionStorage.setItem("rep1_is_pinned_after_reload", String(isPinnedAfterReload));
                sessionStorage.setItem("rep1_note_after_reload", noteAfterReload);

                // Navigate to Report 2
                window.location.href = "/report2.html";
                return;
            }
        } catch (e) {
            var errMarker = document.createElement("div");
            errMarker.id = "browser-probe-result";
            errMarker.textContent = JSON.stringify({ error: e.message, stack: e.stack });
            document.body.appendChild(errMarker);
        }
    }, 250);
});
</script>
"""
        rep2_probe = """
<script>
window.addEventListener("load", function() {
    setTimeout(function() {
        try {
            // Verify Report 2 has 0 pins and no review state from Report 1
            var countBadge = document.getElementById("pinned-count-badge");
            var rep2InitBadge = countBadge ? countBadge.textContent.trim() : "";
            var rep2InitialPinnedCount = window.getPinnedVariantCount();

            // Pin EGFR on Report 2
            var tr0 = $("#variants_table tbody tr[data-variant-key]").first();
            tr0.find("td.dt-control, td.details-control").trigger("click");
            var pinBtn = tr0.next().find(".btn-pin-variant, .pin-toggle-btn");
            pinBtn.trigger("click");
            var noteInp = tr0.next().find(".review-notes-input");
            noteInp.val("EGFR candidate in Report 2");
            noteInp[0].dispatchEvent(new Event("input", { bubbles: true }));

            var rep2BadgeAfterPin = countBadge ? countBadge.textContent.trim() : "";

            // Check localStorage keys
            var rep1StorageKey = "vc_review_report_alpha";
            var rep2StorageKey = window.storageKey; // vc_review_report_beta
            var rep1Raw = window.safeStorageGet(rep1StorageKey);
            var rep2Raw = window.safeStorageGet(rep2StorageKey);

            var result = {
                rep1BadgeStep1: sessionStorage.getItem("rep1_badge_step1"),
                rep1BadgeAfterReload: sessionStorage.getItem("rep1_badge_after_reload"),
                rep1IsPinnedAfterReload: (
                    sessionStorage.getItem("rep1_is_pinned_after_reload") === "true"
                ),
                rep1NoteAfterReload: sessionStorage.getItem("rep1_note_after_reload"),
                rep2InitBadge: rep2InitBadge,
                rep2InitialPinnedCount: rep2InitialPinnedCount,
                rep2BadgeAfterPin: rep2BadgeAfterPin,
                rep1KeyExists: rep1Raw !== null,
                rep2KeyExists: rep2Raw !== null,
                keysDiffer: rep1StorageKey !== rep2StorageKey,
                rep1HasBrca: Boolean(
                    rep1Raw && rep1Raw.indexOf("Persistent review note for BRCA1") !== -1
                ),
                rep2HasEgfr: Boolean(
                    rep2Raw && rep2Raw.indexOf("EGFR candidate in Report 2") !== -1
                ),
                rep2DoesNotHaveBrca: Boolean(rep2Raw && rep2Raw.indexOf("BRCA1") === -1)
            };

            var marker = document.createElement("div");
            marker.id = "browser-probe-result";
            marker.textContent = JSON.stringify(result);
            document.body.appendChild(marker);
        } catch (e) {
            var errMarker = document.createElement("div");
            errMarker.id = "browser-probe-result";
            errMarker.textContent = JSON.stringify({ error: e.message, stack: e.stack });
            document.body.appendChild(errMarker);
        }
    }, 250);
});
</script>
"""
        probed_rep1_html = rep1_html.replace("</body>", rep1_probe + "</body>")
        probed_rep2_html = rep2_html.replace("</body>", rep2_probe + "</body>")

        class MultiReportHandler(http.server.BaseHTTPRequestHandler):
            def do_GET(self):
                self.send_response(200)
                self.send_header("Content-Type", "text/html; charset=utf-8")
                self.end_headers()
                if self.path == "/report2.html":
                    self.wfile.write(probed_rep2_html.encode("utf-8"))
                else:
                    self.wfile.write(probed_rep1_html.encode("utf-8"))

            def log_message(self, *args):
                pass

        server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), MultiReportHandler)
        port = server.server_port
        t = threading.Thread(target=server.serve_forever)
        t.daemon = True
        t.start()

        try:
            res = _dump_dom(chrome_bin, f"http://127.0.0.1:{port}/")
            assert res.returncode == 0, f"Chrome failed: {res.stderr}"
            m = re.search(r'<div id="browser-probe-result">([\s\S]*?)</div>', res.stdout)
            assert m is not None, f"Marker not found: {res.stdout[:500]}"
            import html as html_lib

            data = json.loads(html_lib.unescape(m.group(1)))
            assert "error" not in data, f"Browser error: {data.get('error')}\n{data.get('stack')}"

            # Verify reload persistence on Report 1
            assert data["rep1BadgeStep1"] == "1"
            assert data["rep1BadgeAfterReload"] == "1"
            assert data["rep1IsPinnedAfterReload"] is True
            assert data["rep1NoteAfterReload"] == "Persistent review note for BRCA1 across reload"

            # Verify Report 2 isolation
            assert data["rep2InitBadge"] == "0"
            assert data["rep2InitialPinnedCount"] == 0
            assert data["rep2BadgeAfterPin"] == "1"

            # Verify localStorage separation
            assert data["rep1KeyExists"] is True
            assert data["rep2KeyExists"] is True
            assert data["keysDiffer"] is True
            assert data["rep1HasBrca"] is True
            assert data["rep2HasEgfr"] is True
            assert data["rep2DoesNotHaveBrca"] is True
        finally:
            server.shutdown()
            server.server_close()
