"""Unit tests for HTML report XSS prevention, URL safety, and cell content extraction."""

import json
import re
import subprocess
from pathlib import Path

import pytest


def _load_template_js():
    """Extract helper functions from index.html template for Node execution."""
    template_path = (
        Path(__file__).parent.parent.parent / "variantcentrifuge" / "templates" / "index.html"
    )
    content = template_path.read_text(encoding="utf-8")

    helpers = ""
    for fn_name in ["escapeHtml", "isSafeUrl", "extractCellContent", "formatChildRow"]:
        m = re.search(rf"(function {fn_name}\(.*?\n    \}})", content, re.DOTALL)
        assert m is not None, f"Function {fn_name} not found in index.html"
        helpers += m.group(1) + "\n\n"

    return helpers


def _run_node_script(js_code: str):
    """Run a JS script in node and return stdout."""
    helpers = _load_template_js()
    full_script = f"{helpers}\n{js_code}"
    result = subprocess.run(
        ["node", "-e", full_script],
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout.strip()


@pytest.mark.unit
class TestIsSafeUrl:
    """Test URL safety validator in JavaScript."""

    def test_safe_urls_accepted(self):
        script = """
        const cases = [
            'https://example.org/path?a=1&b=2',
            'http://example.org',
            'igv_reports/chr1:100-200.html',
            'reports/sample1.html',
            'relative_file.html'
        ];
        for (const u of cases) {
            if (!isSafeUrl(u)) {
                console.error('Expected safe: ' + u);
                process.exit(1);
            }
        }
        console.log('OK');
        """
        assert _run_node_script(script) == "OK"

    def test_unsafe_urls_rejected(self):
        script = """
        const cases = [
            'javascript:alert(1)',
            'JAVASCRIPT:alert(1)',
            'data:text/html,<script>alert(1)</script>',
            'vbscript:msgbox(1)',
            '//evil.com/phish',
            'https://example.org/" onclick="alert(1)',
            'https://example.org/<script>',
            'https://example.org/path with spaces',
            'https://example.org/\\\\evil.com',
            'custom:protocol',
            '',
            null,
            undefined
        ];
        for (const u of cases) {
            if (isSafeUrl(u)) {
                console.error('Expected unsafe: ' + u);
                process.exit(1);
            }
        }
        console.log('OK');
        """
        assert _run_node_script(script) == "OK"


@pytest.mark.unit
class TestExtractCellContent:
    """Test cell content extraction from plain text and DOM markup."""

    def test_plain_text(self):
        script = """
        const res = extractCellContent('HIGH');
        if (res.text !== 'HIGH' || res.href !== null || res.links.length !== 0) {
            console.error('Mismatch: ' + JSON.stringify(res));
            process.exit(1);
        }
        console.log('OK');
        """
        assert _run_node_script(script) == "OK"

    def test_truncated_cell_markup(self):
        script = """
        const raw = (
            '<span class="truncated-cell" tabindex="0" ' +
            'data-tippy-content="c.123A>G">c.123A&gt;G</span>'
        );
        const res = extractCellContent(raw);
        if (res.text !== 'c.123A>G' || res.href !== null) {
            console.error('Mismatch: ' + JSON.stringify(res));
            process.exit(1);
        }
        console.log('OK');
        """
        assert _run_node_script(script) == "OK"

    def test_single_anchor(self):
        script = """
        const raw = (
            '<a href="https://example.org" target="_blank">' +
            '<svg class="icon"></svg>' +
            '<span class="sr-only">gnomAD (opens in new tab)</span>gnomAD</a>'
        );
        const res = extractCellContent(raw);
        if (res.text !== 'gnomAD' || res.href !== 'https://example.org' || res.links.length !== 1) {
            console.error('Mismatch: ' + JSON.stringify(res));
            process.exit(1);
        }
        console.log('OK');
        """
        assert _run_node_script(script) == "OK"

    def test_multiple_anchors(self):
        script = """
        const raw = '<a href="igv_reports/s1.html">s1</a> <a href="igv_reports/s2.html">s2</a>';
        const res = extractCellContent(raw);
        if (
            res.links.length !== 2 ||
            res.links[0].href !== 'igv_reports/s1.html' ||
            res.links[1].href !== 'igv_reports/s2.html'
        ) {
            console.error('Mismatch: ' + JSON.stringify(res));
            process.exit(1);
        }
        console.log('OK');
        """
        assert _run_node_script(script) == "OK"


@pytest.mark.unit
class TestFormatChildRow:
    """Test formatChildRow detail panel generation."""

    def test_ordinary_fields_not_rendered_as_links(self):
        """Verify ordinary strings like BRCA1, 100, HIGH are not turned into links."""
        script = """
        const columns = [
            { original_name: 'GENE', display_name: 'GENE' },
            { original_name: 'POS', display_name: 'POS' },
            { original_name: 'IMPACT', display_name: 'IMPACT' }
        ];
        const rowData = [null, 'BRCA1', '100', 'HIGH'];
        const res = formatChildRow(rowData, columns);

        if (res.indexOf('<a') !== -1) {
            console.error('Unexpected anchor in formatChildRow: ' + res);
            process.exit(1);
        }
        if (res.indexOf('<dd>BRCA1</dd>') === -1) {
            console.error('Missing expected GENE text: ' + res);
            process.exit(1);
        }
        if (res.indexOf('<dd>100</dd>') === -1) {
            console.error('Missing expected POS text: ' + res);
            process.exit(1);
        }
        if (res.indexOf('<dd>HIGH</dd>') === -1) {
            console.error('Missing expected IMPACT text: ' + res);
            process.exit(1);
        }
        console.log('OK');
        """
        assert _run_node_script(script) == "OK"

    def test_truncated_hgvs_no_double_escaping(self):
        """Verify truncated cell markup is not double-escaped in detail panel."""
        script = """
        const columns = [
            { original_name: 'HGVS_C', display_name: 'HGVS C' }
        ];
        const raw = (
            '<span class="truncated-cell" data-tippy-content="c.123A>G">c.123A&gt;G</span>'
        );
        const rowData = [null, raw];
        const res = formatChildRow(rowData, columns);

        if (res.indexOf('&lt;span') !== -1 || res.indexOf('<span') !== -1) {
            console.error('Double-escaped or leaked span tag: ' + res);
            process.exit(1);
        }
        if (res.indexOf('<dd>c.123A&gt;G</dd>') === -1) {
            console.error('Expected clean escaped entity c.123A&gt;G: ' + res);
            process.exit(1);
        }
        console.log('OK');
        """
        assert _run_node_script(script) == "OK"

    def test_link_columns_render_proper_anchors(self):
        """Verify link columns and igv_links render safe clickable anchors."""
        script = """
        const columns = [
            {
                original_name: 'gnomAD_link',
                display_name: 'gnomAD',
                is_standard_link_column: true,
                link_display_text: 'gnomAD'
            },
            { original_name: 'igv_links', display_name: 'IGV Reports', is_igv_link_column: true }
        ];
        const rowData = [
            null,
            '<a href="https://example.org" target="_blank">gnomAD</a>',
            '<a href="igv_reports/chr1:100-200.html">sample1</a>'
        ];
        const res = formatChildRow(rowData, columns);

        if (res.indexOf('href="https://example.org"') === -1) {
            console.error('Missing standard link href: ' + res);
            process.exit(1);
        }
        if (res.indexOf('href="igv_reports/chr1:100-200.html"') === -1) {
            console.error('Missing igv link href: ' + res);
            process.exit(1);
        }
        console.log('OK');
        """
        assert _run_node_script(script) == "OK"

    def test_xss_payloads_remain_inert(self):
        """Verify javascript: and injection payloads are neutralised in formatChildRow."""
        script = """
        const columns = [
            { original_name: 'bad_link', display_name: 'Bad Link', is_standard_link_column: true },
            { original_name: 'script_inj', display_name: 'Script Injection' },
            { original_name: 'encoded_inj', display_name: 'Encoded Injection' }
        ];
        const rowData = [
            null,
            'javascript:alert(1)',
            '<script>alert("pwned")</script>',
            '&lt;script&gt;alert("pwned")&lt;/script&gt;'
        ];
        const res = formatChildRow(rowData, columns);

        if (res.indexOf('href="javascript:') !== -1) {
            console.error('javascript: rendered in href: ' + res);
            process.exit(1);
        }
        if (res.indexOf('<script>') !== -1) {
            console.error('Raw <script> tag unescaped: ' + res);
            process.exit(1);
        }
        if (res.indexOf('&lt;script&gt;alert(&quot;pwned&quot;)&lt;/script&gt;') === -1) {
            console.error('Expected escaped script tags for encoded_inj: ' + res);
            process.exit(1);
        }
        console.log('OK');
        """
        assert _run_node_script(script) == "OK"


def _run_column_defs_script(column_data: list, test_body: str):
    """Run a script that initializes columnDefs in node and executes assertions."""
    helpers = _load_template_js()
    template_path = (
        Path(__file__).parent.parent.parent / "variantcentrifuge" / "templates" / "index.html"
    )
    content = template_path.read_text(encoding="utf-8")
    m = re.search(r"(columnData\.forEach\(function\(col, index\) \{[\s\S]*?\n    \}\);)", content)
    assert m is not None, "columnData.forEach loop not found in index.html"
    loop_code = m.group(1)

    full_script = f"""
{helpers}
const dtColumnDefs = [];
const columnData = {json.dumps(column_data)};
{loop_code}
{test_body}
"""
    result = subprocess.run(
        ["node", "-e", full_script],
        capture_output=True,
        text=True,
        check=True,
    )
    return result.stdout.strip()


@pytest.mark.unit
class TestColumnRenderers:
    """Test columnDefs render functions for badges, HGVS, and scores."""

    def test_impact_badge(self):
        script = """
        const defObj = dtColumnDefs.find(d => d.render);
        const render = defObj.render;
        const res1 = render('HIGH', 'display', {});
        const res2 = render('<span class="badge">HIGH</span>', 'display', {});
        const res3 = render('HIGH', 'sort', {});

        if (res1.indexOf('background-color:#dc3545') === -1 || res1.indexOf('>HIGH<') === -1) {
            console.error('Invalid res1: ' + res1);
            process.exit(1);
        }
        if (res2 !== res1) {
            console.error('res2 should match res1 without double-escaping: ' + res2);
            process.exit(1);
        }
        if (res3 !== 'HIGH') {
            console.error('Sort should return clean string: ' + res3);
            process.exit(1);
        }
        console.log('OK');
        """
        assert _run_column_defs_script([{"original_name": "IMPACT"}], script) == "OK"

    def test_clinvar_badge(self):
        script = """
        const defObj = dtColumnDefs.find(d => d.render);
        const render = defObj.render;
        const res1 = render('Pathogenic', 'display', {});
        const res2 = render('likely_pathogenic', 'display', {});
        const res3 = render('conflicting', 'display', {});

        const hasPathogenic = (
            res1.indexOf('background-color:#dc3545') !== -1 &&
            res1.indexOf('>Pathogenic<') !== -1
        );
        if (!hasPathogenic) {
            console.error('Invalid pathogenic: ' + res1);
            process.exit(1);
        }
        const hasLikelyPath = (
            res2.indexOf('background-color:#c0392b') !== -1 &&
            res2.indexOf('>likely pathogenic<') !== -1
        );
        if (!hasLikelyPath) {
            console.error('Invalid likely_pathogenic: ' + res2);
            process.exit(1);
        }
        if (res3.indexOf('background-color:#bf5600') === -1) {
            console.error('Invalid conflicting: ' + res3);
            process.exit(1);
        }
        console.log('OK');
        """
        assert _run_column_defs_script([{"original_name": "ClinVar_CLNSIG"}], script) == "OK"

    def test_hgvs_rendering(self):
        script = """
        const defObj = dtColumnDefs.find(d => d.render);
        const render = defObj.render;
        const shortHgvs = 'c.123A>G';
        const longHgvs = 'c.123456789012345678901234567890A>G';

        const resShort = render(shortHgvs, 'display', {});
        if (resShort !== 'c.123A&gt;G') {
            console.error('Short HGVS escaped mismatch: ' + resShort);
            process.exit(1);
        }

        const resLong = render(longHgvs, 'display', {});
        if (resLong.indexOf('class="truncated-cell"') === -1 || resLong.indexOf('...') === -1) {
            console.error('Long HGVS should be truncated: ' + resLong);
            process.exit(1);
        }

        // Test with already wrapped cell
        const wrapped = (
            '<span class="truncated-cell" data-tippy-content="' + shortHgvs + '">c.123A&gt;G</span>'
        );
        const resWrapped = render(wrapped, 'display', {});
        if (resWrapped !== 'c.123A&gt;G') {
            console.error('Wrapped short HGVS should extract clean text: ' + resWrapped);
            process.exit(1);
        }
        console.log('OK');
        """
        assert _run_column_defs_script([{"original_name": "HGVS_C"}], script) == "OK"

    def test_score_rendering(self):
        script = """
        const defObj = dtColumnDefs.find(d => d.render);
        const render = defObj.render;
        const res1 = render('15.2000', 'display', {});
        const rawScore = (
            '<span class="truncated-cell" data-tippy-content="0.0010">0.0010</span>'
        );
        const res2 = render(rawScore, 'display', {});

        if (res1 !== '15.200') {
            console.error('Score formatting mismatch: ' + res1);
            process.exit(1);
        }
        if (res2 !== '0.001') {
            console.error('Score from wrapped cell mismatch: ' + res2);
            process.exit(1);
        }
        console.log('OK');
        """
        assert _run_column_defs_script([{"original_name": "CADD_phred"}], script) == "OK"
