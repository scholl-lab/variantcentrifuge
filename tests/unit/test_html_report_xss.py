"""Unit tests for HTML report XSS prevention, URL safety, and cell content extraction."""

import http.server
import json
import re
import shutil
import subprocess
import threading
from pathlib import Path

import pytest

DOMPARSER_POLYFILL_JS = r"""
function _decodeEntities(s) {
    return s.replace(/&amp;/g, '&')
            .replace(/&lt;/g, '<')
            .replace(/&gt;/g, '>')
            .replace(/&quot;/g, '"')
            .replace(/&#39;/g, "'")
            .replace(/&#x2F;/g, '/');
}

class DOMNode {
    constructor(nodeType, tagName = '', text = '') {
        this.nodeType = nodeType;
        this.tagName = tagName.toUpperCase();
        this.text = text;
        this.attributes = {};
        this.children = [];
        this.parentNode = null;
    }
    getAttribute(name) {
        const val = this.attributes[name.toLowerCase()];
        return val !== undefined ? val : null;
    }
    setAttribute(name, val) {
        this.attributes[name.toLowerCase()] = String(val);
    }
    appendChild(child) {
        child.parentNode = this;
        this.children.push(child);
    }
    removeChild(child) {
        const idx = this.children.indexOf(child);
        if (idx !== -1) {
            this.children.splice(idx, 1);
            child.parentNode = null;
        }
        return child;
    }
    get textContent() {
        if (this.nodeType === 3) return this.text;
        return this.children.map(c => c.textContent).join('');
    }
    get innerText() {
        return this.textContent;
    }
    cloneNode(deep) {
        const clone = new DOMNode(this.nodeType, this.tagName, this.text);
        clone.attributes = { ...this.attributes };
        if (deep) {
            for (const child of this.children) {
                clone.appendChild(child.cloneNode(true));
            }
        }
        return clone;
    }
    querySelectorAll(selector) {
        const results = [];
        const selectors = selector.split(',').map(s => s.trim().toLowerCase());
        function match(node) {
            if (node.nodeType !== 1) return false;
            for (const sel of selectors) {
                if (sel === 'a' && node.tagName === 'A') return true;
                if (sel === 'svg' && node.tagName === 'SVG') return true;
                const cls = (node.getAttribute('class') || '').split(/\s+/);
                if (sel === '.sr-only' && cls.includes('sr-only')) return true;
                if (sel.startsWith('[') && sel.endsWith(']')) {
                    const attr = sel.slice(1, -1);
                    if (node.getAttribute(attr) !== null) return true;
                }
            }
            return false;
        }
        function traverse(node) {
            for (const child of node.children) {
                if (match(child)) results.push(child);
                traverse(child);
            }
        }
        traverse(this);
        return results;
    }
    querySelector(selector) {
        const all = this.querySelectorAll(selector);
        return all.length > 0 ? all[0] : null;
    }
}

function _parseHtmlToDom(html) {
    const root = new DOMNode(1, 'BODY');
    const stack = [root];
    const tagRegex =
        /<!--[\s\S]*?-->|<(\/)?([a-zA-Z0-9\-]+)((?:\s+[a-zA-Z0-9_\-]+(?:\s*=\s*(?:"[^"]*"|'[^']*'|[^\s>]+))?)*)\s*(\/)?>|([^<]+)/g;
    let match;
    while ((match = tagRegex.exec(html)) !== null) {
        if (match[0].startsWith('<!--')) continue;
        if (match[5]) {
            const text = _decodeEntities(match[5]);
            stack[stack.length - 1].appendChild(new DOMNode(3, '', text));
        } else if (match[2]) {
            const isClosing = Boolean(match[1]);
            const tagName = match[2];
            const attrStr = match[3] || '';
            const isSelfClosing =
                Boolean(match[4]) || ['img', 'br', 'hr', 'input'].includes(tagName.toLowerCase());
            if (isClosing) {
                for (let i = stack.length - 1; i > 0; i--) {
                    if (stack[i].tagName === tagName.toUpperCase()) {
                        stack.splice(i);
                        break;
                    }
                }
            } else {
                const el = new DOMNode(1, tagName);
                const attrRegex =
                    /([a-zA-Z0-9_\-]+)(?:\s*=\s*(?:"([^"]*)"|'([^']*)'|([^\s>]+)))?/g;
                let attrMatch;
                while ((attrMatch = attrRegex.exec(attrStr)) !== null) {
                    const name = attrMatch[1];
                    const val =
                        attrMatch[2] !== undefined
                            ? attrMatch[2]
                            : attrMatch[3] !== undefined
                            ? attrMatch[3]
                            : attrMatch[4] !== undefined
                            ? attrMatch[4]
                            : '';
                    el.setAttribute(name, _decodeEntities(val));
                }
                stack[stack.length - 1].appendChild(el);
                if (!isSelfClosing) {
                    stack.push(el);
                }
            }
        }
    }
    return root;
}

class DOMParser {
    parseFromString(str, mimeType) {
        return { body: _parseHtmlToDom(str) };
    }
}
globalThis.DOMParser = DOMParser;
"""


def _load_template_js():
    """Extract helper functions from index.html template for Node execution."""
    template_path = (
        Path(__file__).parent.parent.parent / "variantcentrifuge" / "templates" / "index.html"
    )
    content = template_path.read_text(encoding="utf-8")

    helpers = ""
    for fn_name in [
        "escapeHtml",
        "encodeSymbolicAlleles",
        "isSafeUrl",
        "extractCellContent",
        "formatChildRow",
    ]:
        m = re.search(rf"(function {fn_name}\(.*?\n    \}})", content, re.DOTALL)
        assert m is not None, f"Function {fn_name} not found in index.html"
        helpers += m.group(1) + "\n\n"

    return helpers


def _run_node_script(js_code: str):
    """Run a JS script in node with inert DOMParser available and return stdout."""
    helpers = _load_template_js()
    full_script = f"{DOMPARSER_POLYFILL_JS}\n{helpers}\n{js_code}"
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

    def test_inheritance_pattern_rendering(self):
        script = """
        const defObj = dtColumnDefs.find(d => d.render);
        const render = defObj.render;
        const resNoneDisp = render('none', 'display', {});
        const resRefDisp = render('reference', 'display', {});
        const resDenovoDisp = render('de_novo', 'display', {});
        const resNoneFilter = render('none', 'filter', {});
        const resRefFilter = render('reference', 'filter', {});
        const resDenovoFilter = render('de_novo', 'filter', {});

        if (resNoneDisp.indexOf('Unknown') === -1) {
            console.error('none display should be Unknown badge: ' + resNoneDisp);
            process.exit(1);
        }
        if (resRefDisp.indexOf('Unknown') === -1) {
            console.error('reference display should be Unknown badge: ' + resRefDisp);
            process.exit(1);
        }
        if (resDenovoDisp.indexOf('de novo') === -1) {
            console.error('de_novo display should be de novo badge: ' + resDenovoDisp);
            process.exit(1);
        }
        if (resNoneFilter !== 'none') {
            console.error('none filter should be "none": ' + resNoneFilter);
            process.exit(1);
        }
        if (resRefFilter !== 'reference') {
            console.error('reference filter should be "reference": ' + resRefFilter);
            process.exit(1);
        }
        if (resDenovoFilter !== 'de_novo') {
            console.error('de_novo filter should be "de_novo": ' + resDenovoFilter);
            process.exit(1);
        }
        console.log('OK');
        """
        assert _run_column_defs_script([{"original_name": "Inheritance_Pattern"}], script) == "OK"

    def test_spliceai_link_column_not_overwritten_by_score_renderer(self):
        script = """
        const defObj = dtColumnDefs.find(d => d.targets === 1 && d.render);
        // For standard link column, no numeric score renderer should overwrite the link
        if (defObj) {
            console.error('SpliceAI link column should NOT have a numeric score renderer');
            process.exit(1);
        }
        console.log('OK');
        """
        col_data = [
            {
                "original_name": "SpliceAI",
                "is_standard_link_column": True,
                "is_igv_link_column": False,
            }
        ]
        assert _run_column_defs_script(col_data, script) == "OK"

    @pytest.mark.parametrize(
        "col_name",
        ["IMPACT", "ClinVar_CLNSIG", "HGVS_C", "VAR_ID", "CADD_phred"],
    )
    def test_renderers_escape_filter_mode(self, col_name: str):
        """Verify renderers escape HTML when type='filter' (#130)."""
        xss_payload = "&lt;img src=x onerror=alert(1)&gt;"
        script = f"""
        const defObj = dtColumnDefs.find(d => d.render);
        const render = defObj.render;
        const res = render({json.dumps(xss_payload)}, 'filter', {{}});
        if (res.indexOf('<img') !== -1) {{
            console.error('Filter mode returned unescaped HTML: ' + res);
            process.exit(1);
        }}
        if (res.indexOf('&lt;img') === -1) {{
            console.error('Filter mode did not escape payload: ' + res);
            process.exit(1);
        }}
        console.log('OK');
        """
        assert _run_column_defs_script([{"original_name": col_name}], script) == "OK"

    def test_categorical_filter_option_extraction_and_matching_with_ampersand(self):
        """Verify categorical filter decodes &amp; for dropdown options and filterFn (#130)."""
        script = """
        function isMissing(val) {
            return (
                val === null ||
                val === undefined ||
                val === '' ||
                val === '-' ||
                val === 'NA' ||
                val === 'N/A'
            );
        }
        const rawCatValues = [
            'Pathogenic &amp; uncertain',
            'Benign',
            'Pathogenic &amp; uncertain'
        ];
        const uniqueValues = [];
        const seen = {};
        rawCatValues.forEach(function(v) {
            if (!isMissing(v)) {
                const str = String(extractCellContent(v).text || '').trim();
                if (str && !seen[str]) {
                    seen[str] = true;
                    uniqueValues.push(str);
                }
            }
        });
        uniqueValues.sort();

        if (uniqueValues.length !== 2 || uniqueValues[1] !== 'Pathogenic & uncertain') {
            console.error('Dropdown uniqueValues mismatch: ' + JSON.stringify(uniqueValues));
            process.exit(1);
        }

        const selectedValue = 'Pathogenic & uncertain';
        const filterFn = function(settings, data) {
            if (settings.nTable.id !== 'variants_table') return true;
            const cellValue = data[0];
            if (isMissing(cellValue)) return false;
            return String(cellValue).trim() === selectedValue;
        };

        const settings = { nTable: { id: 'variants_table' } };
        if (!filterFn(settings, ['Pathogenic & uncertain'])) {
            console.error('filterFn failed to match decoded text');
            process.exit(1);
        }
        if (filterFn(settings, ['Benign'])) {
            console.error('filterFn incorrectly matched Benign');
            process.exit(1);
        }

        const lowerSearch = 'pathogenic & uncertain';
        const textFilterFn = function(settings, data) {
            if (settings.nTable.id !== 'variants_table') return true;
            const cellValue = data[0];
            if (isMissing(cellValue)) return false;
            return String(cellValue).toLowerCase().indexOf(lowerSearch) !== -1;
        };
        if (!textFilterFn(settings, ['Pathogenic & uncertain'])) {
            console.error('textFilterFn failed to match plain text');
            process.exit(1);
        }
        console.log('OK');
        """
        assert _run_node_script(script) == "OK"

    def test_extract_cell_content_uses_dom_parser_and_filter_fn_plain_text(self):
        """Verify extractCellContent uses DOMParser and filterFn compares plain text (#130)."""
        template_path = (
            Path(__file__).parent.parent.parent / "variantcentrifuge" / "templates" / "index.html"
        )
        content = template_path.read_text(encoding="utf-8")

        # 1. extractCellContent should use DOMParser, never document.createElement('div').innerHTML
        m_fn = re.search(r"function extractCellContent\(raw\) \{([\s\S]*?)\n    \}", content)
        assert m_fn is not None
        fn_body = m_fn.group(1)
        assert "DOMParser" in fn_body
        assert "document.createElement('div')" not in fn_body
        assert "div.innerHTML" not in fn_body

        # 2. filterFn should compare plain cellValue, not call extractCellContent
        m_filters = re.search(
            r"else if \(filterType === 'categorical'\) \{([\s\S]*?)\n    \} // end if filterRow",
            content,
        )
        assert m_filters is not None
        filters_body = m_filters.group(1)
        assert "String(cellValue).trim() === selectedValue" in filters_body
        assert "extractCellContent(cellValue)" not in filters_body

    def test_dom_parser_pipeline_and_xss_mitigations(self):
        """Exercise DOMParser extractCellContent, rendering, sorting, and filtering (#130)."""
        script = """
        // 1. extractCellContent executes DOMParser branch
        const payload = (
            '<a href="javascript:alert(1)"><svg><circle></circle></svg>' +
            '<span class="sr-only">hide</span>safe_label</a>'
        );
        const extracted = extractCellContent(payload);
        if (extracted.text !== 'safe_label' || extracted.href !== 'javascript:alert(1)') {
            console.error('DOMParser extractCellContent failed: ' + JSON.stringify(extracted));
            process.exit(1);
        }
        if (isSafeUrl(extracted.href)) {
            console.error('javascript: URL was unexpectedly marked safe');
            process.exit(1);
        }

        // 2. Column rendering: IMPACT badge escaping
        const impactRender = function(data, type, row) {
            if (data) {
                var cell = extractCellContent(data);
                if (type === 'display') {
                    var colors = { 'HIGH': '#dc3545', 'MODERATE': '#c05000' };
                    var color = colors[cell.text] || '#6c757d';
                    return '<span class="badge" style="background-color:' + color + ';">' +
                        escapeHtml(cell.text) + '</span>';
                }
                return cell.text ? escapeHtml(String(cell.text)) : '';
            }
            return data != null ? escapeHtml(String(data)) : '';
        };

        const impactPayload = '<script>alert("impact")</script>';
        const renderedImpact = impactRender(impactPayload, 'display', []);
        if (
            renderedImpact.indexOf('<script>') !== -1 ||
            (renderedImpact.indexOf('&lt;script&gt;') === -1 &&
             renderedImpact.indexOf('alert(&quot;impact&quot;)') === -1)
        ) {
            console.error('Impact render did not escape: ' + renderedImpact);
            process.exit(1);
        }

        // 3. Sorting (type === 'sort'): Returns plain clean text, not HTML markup
        const sortedImpact = impactRender('<span class="badge">HIGH</span>', 'sort', []);
        if (sortedImpact !== 'HIGH') {
            console.error('Sort render failed to extract plain text: ' + sortedImpact);
            process.exit(1);
        }

        // 4. Tooltip attribute escaping: HGVS column > 30 chars
        const hgvsRender = function(data, type, row) {
            if (data) {
                var cell = extractCellContent(data);
                var str = String(cell.text || '');
                if (type === 'display') {
                    if (str.length > 30) {
                        var truncated = str.substring(0, 30) + '...';
                        return '<span class="truncated-cell" tabindex="0" data-tippy-content="' +
                            escapeHtml(str) + '">' + escapeHtml(truncated) + '</span>';
                    }
                    return escapeHtml(str);
                }
                return escapeHtml(str);
            }
            return data != null ? escapeHtml(String(data)) : '';
        };

        const hgvsPayload = (
            'c.123A>G" onfocus="alert(1)" data-breakout="excessively_long_variant_hgvs_string'
        );
        const renderedHgvs = hgvsRender(hgvsPayload, 'display', []);
        if (
            renderedHgvs.indexOf('" onfocus=') !== -1 ||
            renderedHgvs.indexOf('&quot; onfocus=&quot;') === -1
        ) {
            console.error('Tooltip attribute not properly escaped: ' + renderedHgvs);
            process.exit(1);
        }
        const sortedHgvs = hgvsRender(hgvsPayload, 'sort', []);
        if (sortedHgvs.indexOf('<span') !== -1) {
            console.error('Sort mode returned markup for HGVS: ' + sortedHgvs);
            process.exit(1);
        }

        // 5. Filtering (filterFn): compares plain cellValue without DOM re-parsing
        const filterFn = function(settings, data) {
            if (settings.nTable.id !== 'variants_table') return true;
            var cellValue = data[0];
            return String(cellValue).trim() === 'expected_value';
        };
        const settings = { nTable: { id: 'variants_table' } };
        if (filterFn(settings, [hgvsPayload]) !== false) {
            console.error('filterFn matched unexpected payload');
            process.exit(1);
        }

        // 6. Child-row details rendering: formatChildRow against XSS
        const columns = [
            { original_name: 'inj_col', display_name: 'Injection' },
            { original_name: 'bad_link', display_name: 'Bad Link', is_standard_link_column: true }
        ];
        const rowData = [
            null,
            '<img src=x onerror=alert(1)>',
            '<a href="javascript:alert(2)">Click</a>'
        ];
        const childHtml = formatChildRow(rowData, columns);
        if (childHtml.indexOf('<img') !== -1 || childHtml.indexOf('onerror=') !== -1) {
            console.error('Child row leaked unescaped img/onerror tag: ' + childHtml);
            process.exit(1);
        }
        if (childHtml.indexOf('href="javascript:') !== -1) {
            console.error('Child row rendered dangerous javascript: href: ' + childHtml);
            process.exit(1);
        }

        console.log('OK');
        """
        assert _run_node_script(script) == "OK"

    def test_full_rendered_main_script_initialization_and_interactions(self, tmp_path: Path):
        """Execute full un-sliced main <script> of rendered report in Node (#130)."""
        from variantcentrifuge.generate_html_report import generate_html_report

        variants_data = [
            {
                "CHROM": "chr1",
                "POS": 100,
                "REF": "A",
                "ALT": "T",
                "IMPACT": '<script>alert("impact")</script>',
                "ClinVar_CLNSIG": "Pathogenic & uncertain",
                "Inheritance_Pattern": "de_novo",
                "GENE": "BRCA1",
                "CADD_phred": "25.4",
                "gnomAD_AF": "0.00005",
                "CustomLink": '<a href="javascript:alert(2)">Click</a>',
            },
            {
                "CHROM": "chr1",
                "POS": 200,
                "REF": "C",
                "ALT": "G",
                "IMPACT": "MODERATE",
                "ClinVar_CLNSIG": "Benign",
                "Inheritance_Pattern": "autosomal_recessive",
                "GENE": "TP53",
                "CADD_phred": "12.0",
                "gnomAD_AF": "0.01",
                "CustomLink": '<a href="https://example.org/safe">Safe</a>',
            },
        ]
        summary = {
            "num_variants": 2,
            "num_samples": 1,
            "impact_distribution": {"HIGH": 1, "MODERATE": 1},
            "inheritance_distribution": {"de_novo": 1, "autosomal_recessive": 1},
            "top_genes": [{"gene": "BRCA1", "count": 1}, {"gene": "TP53", "count": 1}],
        }
        cfg = {"reference": "GRCh38", "filters": "PASS"}
        v_file = tmp_path / "variants.json"
        s_file = tmp_path / "summary.json"
        v_file.write_text(json.dumps(variants_data), encoding="utf-8")
        s_file.write_text(json.dumps(summary), encoding="utf-8")
        generate_html_report(str(v_file), str(s_file), str(tmp_path), cfg)
        report_html = (tmp_path / "index.html").read_text(encoding="utf-8")

        all_scripts = re.findall(r"<script>([\s\S]*?)</script>", report_html)
        assert len(all_scripts) == 11, f"Expected 11 <script> blocks, got {len(all_scripts)}"
        bundled_scripts = all_scripts[:-1]
        main_script = all_scripts[-1]

        # 1. Assert that corrupting a bundled script causes Node execution to fail
        corrupted_bundled = list(bundled_scripts)
        corrupted_bundled[2] = "throw new Error('BUNDLED_STARTUP_BOMB');"
        bomb_code = f"""
delete globalThis.module;
delete globalThis.exports;
delete globalThis.require;
{corrupted_bundled[2]}
"""
        res_bomb = subprocess.run(["node"], input=bomb_code, capture_output=True, text=True)
        assert res_bomb.returncode != 0, "Corrupted bundled script should fail in Node"

        # 2. Assert that running all 10 bundled scripts in Node defines window.jQuery
        # and window.jQuery.fn.dataTable
        bundled_env_setup = f"""
{DOMPARSER_POLYFILL_JS}

// Disable Node CommonJS so UMD modules attach to window/globalThis
delete globalThis.module;
delete globalThis.exports;
delete globalThis.require;

function createMockElement(tag, id = '') {{
    const el = {{
        nodeType: 1,
        nodeName: tag.toUpperCase(),
        tagName: tag.toUpperCase(),
        id: id,
        checked: true,
        lastChild: {{ checked: true, nodeType: 1 }},
        firstChild: {{ checked: true, nodeType: 1 }},
        childNodes: [],
        children: [],
        style: {{ display: '', setProperty: () => {{}}, getPropertyValue: () => '' }},
        classList: {{
            _classes: new Set(),
            add(c) {{ this._classes.add(c); }},
            remove(c) {{ this._classes.delete(c); }},
            contains(c) {{ return this._classes.has(c); }},
            toggle(c) {{ if (this.contains(c)) this.remove(c); else this.add(c); }}
        }},
        attributes: {{}},
        setAttribute(k, v) {{ this.attributes[k] = String(v); }},
        getAttribute(k) {{ return this.attributes[k] || null; }},
        appendChild(c) {{
            if (c) {{
                c.parentNode = this;
                this.childNodes.push(c);
                this.children.push(c);
            }}
            return c;
        }},
        removeChild(c) {{
            const idx = this.childNodes.indexOf(c);
            if (idx !== -1) this.childNodes.splice(idx, 1);
            return c;
        }},
        replaceChild(newC, oldC) {{ return newC; }},
        insertBefore(c, ref) {{ return this.appendChild(c); }},
        addEventListener(ev, fn) {{}},
        click() {{}},
        closest(sel) {{ return this; }},
        getContext(type) {{ return {{}}; }},
        querySelector(sel) {{ return createMockElement('div'); }},
        querySelectorAll(sel) {{ return []; }},
        getElementsByTagName(tag) {{ return [createMockElement('td')]; }},
        getBoundingClientRect() {{
            return {{ top: 0, left: 0, right: 0, bottom: 0, width: 0, height: 0 }};
        }},
        cloneNode(deep) {{ return createMockElement(this.tagName); }},
        parentNode: null,
        textContent: '',
        innerHTML: ''
    }};
    el.parentNode = {{
        replaceChild: () => {{}},
        appendChild: () => {{}},
        removeChild: () => {{}}
    }};
    return el;
}}

const elementsById = new Map();
function getOrCreateElement(id) {{
    if (!elementsById.has(id)) {{
        elementsById.set(id, createMockElement('div', id));
    }}
    return elementsById.get(id);
}}

['variants_table', 'table-skeleton', 'toolbar-slot-length', 'toolbar-slot-search',
 'toolbar-slot-buttons', 'active-filters-strip', 'filter-row-top', 'filter-row-bottom',
 'impact_chart', 'variant_type_chart', 'chromosome_chart',
 'af_histogram_chart', 'inheritance_chart', 'variants-table-wrapper',
 'viz-collapse-toggle'].forEach(getOrCreateElement);

globalThis.document = {{
    nodeType: 9,
    documentElement: createMockElement('html'),
    head: createMockElement('head'),
    body: createMockElement('body'),
    createElement: (tag) => createMockElement(tag),
    createElementNS: (ns, t) => createMockElement(t),
    createDocumentFragment: () => createMockElement('#fragment'),
    createTextNode: (t) => ({{ nodeValue: String(t), nodeType: 3 }}),
    createComment: (t) => ({{ nodeValue: String(t), nodeType: 8 }}),
    createEvent: () => ({{ initEvent: () => {{}} }}),
    implementation: {{ createHTMLDocument: () => ({{ body: createMockElement('body') }}) }},
    defaultView: globalThis,
    getElementById: (id) => getOrCreateElement(id),
    querySelector: (sel) => {{
        if (sel.startsWith('#')) return getOrCreateElement(sel.slice(1));
        return createMockElement('div');
    }},
    querySelectorAll: (sel) => [],
    addEventListener: (ev, fn) => {{
        if (ev === 'DOMContentLoaded') {{
            globalThis._domLoadedCallback = fn;
        }}
    }}
}};

globalThis.window = globalThis;
globalThis.window.document = globalThis.document;
globalThis.navigator = {{ userAgent: 'Mozilla/5.0 (X11; Linux x86_64)', platform: 'Linux' }};
globalThis.location = {{
    href: 'http://localhost/',
    protocol: 'http:',
    host: 'localhost',
    hostname: 'localhost'
}};
globalThis.Element = function() {{}};
globalThis.HTMLElement = function() {{}};
globalThis.Node = function() {{}};
globalThis.MutationObserver = function() {{
    this.observe = () => {{}};
    this.disconnect = () => {{}};
}};
globalThis.ResizeObserver = function() {{
    this.observe = () => {{}};
    this.disconnect = () => {{}};
}};
window.addEventListener = (ev, fn) => {{}};
window.print = () => {{}};
window.matchMedia = () => ({{ matches: false, addEventListener: () => {{}} }});

const storage = new Map();
globalThis.localStorage = {{
    getItem: (k) => storage.get(k) || null,
    setItem: (k, v) => storage.set(k, String(v)),
    removeItem: (k) => storage.delete(k),
    clear: () => storage.clear()
}};
"""
        bundled_concat = "\n;\n".join(bundled_scripts)
        verify_bundled_js = """
if (typeof window.jQuery !== 'function') throw new Error('window.jQuery is not a function');
if (typeof (window.jQuery.fn && window.jQuery.fn.dataTable) !== 'function') {
    throw new Error('window.jQuery.fn.dataTable is not a function');
}
console.log('BUNDLED_ASSETS_OK');
"""
        res_bundled = subprocess.run(
            ["node"],
            input=f"{bundled_env_setup}\n{bundled_concat}\n{verify_bundled_js}",
            capture_output=True,
            text=True,
        )
        assert res_bundled.returncode == 0, f"Bundled assets failed: {res_bundled.stderr}"
        assert "BUNDLED_ASSETS_OK" in res_bundled.stdout

        test_rows = [
            [
                None,
                "chr1",
                100,
                "A",
                "T",
                '<script>alert("impact")</script>',
                "Pathogenic & uncertain",
                "de_novo",
                "BRCA1",
                "25.4",
                "0.00005",
                '<a href="javascript:alert(2)">Click</a>',
            ],
            [
                None,
                "chr1",
                200,
                "C",
                "G",
                "MODERATE",
                "Benign",
                "autosomal_recessive",
                "TP53",
                "12.0",
                "0.01",
                '<a href="https://example.org/safe">Safe</a>',
            ],
        ]

        node_harness = f"""
{DOMPARSER_POLYFILL_JS}

const createdSelects = [];

function createMockElement(tag, id = '') {{
    const el = {{
        tagName: tag.toUpperCase(),
        id: id,
        dataset: {{}},
        _listeners: {{}},
        style: {{ display: '', setProperty: () => {{}}, getPropertyValue: () => '' }},
        classList: {{
            _classes: new Set(),
            add(c) {{ this._classes.add(c); }},
            remove(c) {{ this._classes.delete(c); }},
            contains(c) {{ return this._classes.has(c); }},
            toggle(c) {{ if (this.contains(c)) this.remove(c); else this.add(c); }}
        }},
        attributes: {{}},
        setAttribute(k, v) {{ this.attributes[k] = String(v); }},
        getAttribute(k) {{ return this.attributes[k] || null; }},
        appendChild(c) {{
            if (c) {{
                c.parentNode = this;
                if (c.tagName === 'OPTION' && this.options) {{
                    this.options.push(c);
                }}
            }}
            return c;
        }},
        removeChild(c) {{ return c; }},
        replaceChild(newC, oldC) {{ return newC; }},
        insertBefore(c, ref) {{ return this.appendChild(c); }},
        addEventListener(ev, fn) {{
            this._listeners[ev] = fn;
        }},
        click() {{}},
        closest(sel) {{ return this; }},
        getContext(type) {{ return {{}}; }},
        querySelector(sel) {{ return createMockElement('div'); }},
        querySelectorAll(sel) {{ return []; }},
        parentNode: null,
        textContent: '',
        innerHTML: ''
    }};
    if (tag.toLowerCase() === 'select') {{
        el.options = [];
        createdSelects.push(el);
    }}
    el.parentNode = {{
        replaceChild: () => {{}},
        appendChild: () => {{}},
        removeChild: () => {{}}
    }};
    return el;
}}

const elementsById = new Map();
function getOrCreateElement(id) {{
    if (!elementsById.has(id)) {{
        elementsById.set(id, createMockElement('div', id));
    }}
    return elementsById.get(id);
}}

['variants_table', 'table-skeleton', 'toolbar-slot-length', 'toolbar-slot-search',
 'toolbar-slot-buttons', 'active-filters-strip', 'filter-row',
 'impact_chart', 'variant_type_chart', 'chromosome_chart',
 'af_histogram_chart', 'inheritance_chart', 'variants-table-wrapper',
 'viz-collapse-toggle'].forEach(getOrCreateElement);

globalThis.window = globalThis;
window.addEventListener = (ev, fn) => {{}};
window.print = () => {{}};
window.matchMedia = () => ({{ matches: false, addEventListener: () => {{}} }});

const storage = new Map();
globalThis.localStorage = {{
    getItem: (k) => storage.get(k) || null,
    setItem: (k, v) => storage.set(k, String(v)),
    removeItem: (k) => storage.delete(k),
    clear: () => storage.clear()
}};

globalThis.document = {{
    getElementById: (id) => getOrCreateElement(id),
    querySelector: (sel) => {{
        if (sel.startsWith('#')) return getOrCreateElement(sel.slice(1));
        return createMockElement('div');
    }},
    querySelectorAll: (sel) => [],
    createElement: (tag) => createMockElement(tag),
    createTextNode: (t) => createMockElement('#text'),
    body: createMockElement('body'),
    addEventListener: (ev, fn) => {{
        if (ev === 'DOMContentLoaded') {{
            globalThis._domLoadedCallback = fn;
        }}
    }}
}};

globalThis.tippy = function(el, opts) {{
    if (el) el._tippy = {{ destroy: () => {{}} }};
    return el._tippy;
}};

globalThis.noUiSlider = {{
    create: function(el, opts) {{
        el.noUiSlider = {{
            on: function(ev, cb) {{
                if (ev === 'update') cb([opts.range.min, opts.range.max]);
            }},
            set: function(v) {{}},
            get: function() {{ return [opts.range.min, opts.range.max]; }}
        }};
    }}
}};

function MockChart(ctx, config) {{
    this.ctx = ctx;
    this.config = config;
    this.destroy = () => {{}};
    this.update = () => {{}};
}}
globalThis.Chart = MockChart;
globalThis.ChartDataLabels = {{}};

const testRows = {json.dumps(test_rows)};
let registeredTable = null;
let clickHandler = null;

function MockDataTable(sel, options) {{
    registeredTable = this;
    this.options = options;
    this.columnDefs = options.columnDefs || [];

    this.on = function(event, arg2, arg3) {{
        if (event === 'click') {{
            clickHandler = arg3 || arg2;
        }}
    }};

    this.api = function() {{
        return {{
            columns: {{ adjust: () => {{}} }}
        }};
    }};

    this.table = function() {{
        return {{
            node: () => getOrCreateElement('variants_table')
        }};
    }};

    this.columns = {{
        adjust: () => ({{ draw: () => {{}} }})
    }};

    this.column = function(ci) {{
        return {{
            data: () => ({{
                toArray: () => testRows.map(r => r[ci])
            }})
        }};
    }};

    this.row = function(tr) {{
        return {{
            child: Object.assign(function(content) {{
                globalThis._lastChildContent = content;
                return {{ show: () => {{}} }};
            }}, {{
                isShown: () => false,
                hide: () => {{}}
            }}),
            data: () => testRows[0]
        }};
    }};

    this.rows = function(sel) {{
        return {{
            data: () => ({{
                toArray: () => testRows
            }})
        }};
    }};

    this.draw = function() {{}};

    if (options.initComplete) {{
        options.initComplete.call(this);
    }}
    if (options.drawCallback) {{
        options.drawCallback.call(this, {{}});
    }}
}}

MockDataTable.ext = {{ search: [] }};

const jQueryMock = function(sel) {{
    if (typeof sel === 'function') {{
        sel();
        return;
    }}
    return {{
        ready: (cb) => cb(),
        on: () => {{}},
        text: () => {{}}
    }};
}};
jQueryMock.fn = {{
    dataTable: MockDataTable,
    DataTable: MockDataTable
}};
globalThis.$ = jQueryMock;
globalThis.jQuery = jQueryMock;
globalThis.DataTable = MockDataTable;

// Execute the main script!
{main_script}

// Run DOMContentLoaded
if (globalThis._domLoadedCallback) {{
    globalThis._domLoadedCallback();
}}

// Verify DataTable was initialized
if (!registeredTable) {{
    console.error('DataTable was not instantiated');
    process.exit(1);
}}

// 1. Exercise columnDefs renderers for display, sort, and filter on all test rows
for (const colDef of registeredTable.columnDefs) {{
    if (colDef.render) {{
        for (const row of testRows) {{
            const val = row[colDef.targets] !== undefined ? row[colDef.targets] : row[1];
            for (const mode of ['display', 'sort', 'filter']) {{
                const res = colDef.render(val, mode, row);
                if (typeof res === 'string') {{
                    if (res.includes('<script>') || res.includes('onerror=')) {{
                        console.error('XSS payload not escaped by columnDef:', res);
                        process.exit(1);
                    }}
                }}
            }}
        }}
    }}
}}

// 2. Exercise categorical filter selection and registered search callbacks
const clinvarSelect = createdSelects.find(s => s.id === 'filter-ClinVar_CLNSIG');
if (!clinvarSelect) {{
    console.error('ClinVar_CLNSIG select not created');
    process.exit(1);
}}
clinvarSelect.value = 'Pathogenic & uncertain';
if (!clinvarSelect._listeners.change) {{
    console.error('No change listener on clinvarSelect');
    process.exit(1);
}}
clinvarSelect._listeners.change.call(clinvarSelect);

if (MockDataTable.ext.search.length === 0) {{
    console.error('Expected MockDataTable.ext.search to have registered filters');
    process.exit(1);
}}
const registeredFilterFn = MockDataTable.ext.search[MockDataTable.ext.search.length - 1];
const match0 = registeredFilterFn({{ nTable: {{ id: 'variants_table' }} }}, testRows[0]);
const match1 = registeredFilterFn({{ nTable: {{ id: 'variants_table' }} }}, testRows[1]);
if (!match0 || match1) {{
    console.error('Search callback match mismatch: row0=' + match0 + ', row1=' + match1);
    process.exit(1);
}}

// 3. Exercise child-row click handler
if (clickHandler) {{
    const mockEvent = {{
        target: {{
            closest: (sel) => ({{
                classList: {{
                    add: () => {{}},
                    remove: () => {{}}
                }}
            }})
        }}
    }};
    clickHandler(mockEvent);
    if (globalThis._lastChildContent) {{
        const cnt = globalThis._lastChildContent;
        if (cnt.includes('<script>') || cnt.includes('javascript:alert')) {{
            console.error('Child row content leaked dangerous tag or URL:', cnt);
            process.exit(1);
        }}
    }}
}}

console.log('ALL_CHECKS_PASSED');
"""

        res = subprocess.run(["node"], input=node_harness, capture_output=True, text=True)
        assert res.returncode == 0, f"Node execution failed: {res.stderr}"
        assert res.stdout.strip() == "ALL_CHECKS_PASSED"

        # 4. Verify that mutating DOMContentLoaded to throw causes failure
        mut_target = "document.addEventListener('DOMContentLoaded', function() {"
        mut_replacement = (
            "document.addEventListener('DOMContentLoaded', function() { "
            "throw new Error('Startup bomb');"
        )
        mutated_script = main_script.replace(mut_target, mut_replacement)
        mutated_harness = node_harness.replace(main_script, mutated_script)
        res_mutated = subprocess.run(
            ["node"], input=mutated_harness, capture_output=True, text=True
        )
        assert res_mutated.returncode != 0, "Mutated startup should have failed execution"

    def test_generated_report_real_browser_datatables_interactions(self, tmp_path: Path):
        """Serve full rendered HTML report and verify DataTables in headless Chrome (#130)."""
        from variantcentrifuge.generate_html_report import generate_html_report

        chrome_bin = (
            shutil.which("google-chrome")
            or shutil.which("chromium")
            or shutil.which("chromium-browser")
        )
        if not chrome_bin:
            pytest.skip("Chrome/Chromium binary not available")

        variants_data = [
            {
                "CHROM": "chr1",
                "POS": 100,
                "REF": "A",
                "ALT": "T",
                "IMPACT": '<script>alert("impact")</script>',
                "ClinVar_CLNSIG": "Pathogenic & uncertain <img src=x onerror=alert(1)>",
                "Inheritance_Pattern": "de_novo",
                "GENE": "BRCA1",
                "CADD_phred": "25.4",
                "gnomAD_AF": "0.00005",
                "CustomLink": '<a href="javascript:alert(2)">Click</a>',
            },
            {
                "CHROM": "chr1",
                "POS": 200,
                "REF": "C",
                "ALT": "G",
                "IMPACT": "MODERATE",
                "ClinVar_CLNSIG": "Benign",
                "Inheritance_Pattern": "autosomal_recessive",
                "GENE": "TP53",
                "CADD_phred": "12.0",
                "gnomAD_AF": "0.01",
                "CustomLink": '<a href="https://example.org/safe">Safe</a>',
            },
        ]
        summary = {
            "num_variants": 2,
            "num_samples": 1,
            "impact_distribution": {"HIGH": 1, "MODERATE": 1},
            "inheritance_distribution": {"de_novo": 1, "autosomal_recessive": 1},
            "top_genes": [{"gene": "BRCA1", "count": 1}, {"gene": "TP53", "count": 1}],
        }
        cfg = {
            "reference": "GRCh38",
            "filters": "PASS",
            "html_report_truncate_settings": {
                "columns_for_hover_expand": ["ClinVar_CLNSIG"],
                "default_max_width_px": 80,
            },
        }
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
            var initCount = dt.rows().count();

            var headers = [];
            $("#variants_table thead th").each(function() { headers.push($(this).text().trim()); });
            var posIdx = headers.indexOf("POS");
            if (posIdx === -1) posIdx = 2;

            dt.order([posIdx, "desc"]).draw();
            var firstRowPosDesc = $("#variants_table tbody tr")
                .first().find("td").eq(posIdx).text().trim();

            dt.order([posIdx, "asc"]).draw();
            var firstRowPosAsc = $("#variants_table tbody tr")
                .first().find("td").eq(posIdx).text().trim();

            dt.search("BRCA1").draw();
            var filteredCount = dt.rows({ search: "applied" }).count();
            dt.search("").draw();

            var clinvarSelect = $("#filter-ClinVar_CLNSIG");
            var clinvarFilteredCount = -1;
            if (clinvarSelect.length > 0) {
                var targetOpt = clinvarSelect.find("option").filter(function() {
                    return $(this).text().indexOf("Pathogenic & uncertain") !== -1;
                }).val();
                if (targetOpt) {
                    clinvarSelect.val(targetOpt);
                    clinvarSelect[0].dispatchEvent(new Event("change"));
                    clinvarFilteredCount = dt.rows({ search: "applied" }).count();
                    clinvarSelect.val("");
                    clinvarSelect[0].dispatchEvent(new Event("change"));
                }
            }

            var firstControl = $("#variants_table tbody td.dt-control").first();
            firstControl.trigger("click");
            var isShown = dt.row(0).child.isShown();
            var childHtml = dt.row(0).child() ? String($(dt.row(0).child()).html()) : "";
            var hasRawScript = childHtml.indexOf("<script>") !== -1;

            var tippyCells = Array.from(
                document.querySelectorAll(".truncated-cell[data-tippy-content]")
            );
            var tippyInitialized = tippyCells.length > 0 && tippyCells.every(function(c) {
                return Boolean(c._tippy);
            });
            var tippyContent = tippyCells.map(function(c) {
                return c.getAttribute("data-tippy-content") || "";
            });

            var result = {
                initCount: initCount,
                firstRowPosDesc: firstRowPosDesc,
                firstRowPosAsc: firstRowPosAsc,
                filteredCount: filteredCount,
                clinvarFilteredCount: clinvarFilteredCount,
                isShown: isShown,
                hasRawScript: hasRawScript,
                tippyInitialized: tippyInitialized,
                tippyContent: tippyContent
            };
            var marker = document.createElement("div");
            marker.id = "browser-probe-result";
            marker.textContent = JSON.stringify(result);
            document.body.appendChild(marker);
        } catch (e) {
            var errMarker = document.createElement("div");
            errMarker.id = "browser-probe-result";
            errMarker.textContent = JSON.stringify({ error: e.message });
            document.body.appendChild(errMarker);
        }
    }, 200);
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

        server = http.server.HTTPServer(("127.0.0.1", 0), Handler)
        port = server.server_port
        t = threading.Thread(target=server.serve_forever)
        t.daemon = True
        t.start()

        try:
            cmd = [
                chrome_bin,
                "--headless",
                "--disable-gpu",
                "--no-sandbox",
                "--virtual-time-budget=2000",
                "--dump-dom",
                f"http://127.0.0.1:{port}/",
            ]
            res = subprocess.run(cmd, capture_output=True, text=True, timeout=15)
            assert res.returncode == 0, f"Headless Chrome execution failed: {res.stderr}"
            m = re.search(r'<div id="browser-probe-result">([\s\S]*?)</div>', res.stdout)
            assert m is not None, f"Probe marker not found in Chrome DOM output: {res.stdout[:500]}"
            import html as html_lib

            data = json.loads(html_lib.unescape(m.group(1)))
            assert "error" not in data, f"Browser probe reported error: {data.get('error')}"
            assert data["initCount"] == 2
            assert data["firstRowPosDesc"] == "200"
            assert data["firstRowPosAsc"] == "100"
            assert data["filteredCount"] == 1
            assert data["clinvarFilteredCount"] == 1
            assert data["isShown"] is True
            assert data["hasRawScript"] is False
            assert data["tippyInitialized"] is True
            assert "Pathogenic & uncertain <img src=x onerror=alert(1)>" in data["tippyContent"]
        finally:
            server.shutdown()
            server.server_close()

    def test_report_width_configuration_xss_injection_prevented(self, tmp_path: Path):
        """Verify malicious default_max_width_px cannot inject onmouseover handler (#130)."""
        from variantcentrifuge.generate_html_report import generate_html_report

        chrome_bin = (
            shutil.which("google-chrome")
            or shutil.which("chromium")
            or shutil.which("chromium-browser")
        )
        if not chrome_bin:
            pytest.skip("Chrome/Chromium binary not available")

        variants_data = [
            {
                "CHROM": "chr1",
                "POS": 100,
                "REF": "A",
                "ALT": "T",
                "IMPACT": "HIGH",
                "ClinVar_CLNSIG": "Pathogenic",
                "Inheritance_Pattern": "de_novo",
                "GENE": "BRCA1",
            }
        ]
        summary = {
            "num_variants": 1,
            "num_samples": 1,
            "impact_distribution": {"HIGH": 1},
            "inheritance_distribution": {"de_novo": 1},
            "top_genes": [{"gene": "BRCA1", "count": 1}],
        }
        xss_width = '150px;" onmouseover="globalThis.WIDTH_PROBE=1" data-x="'
        cfg = {
            "reference": "GRCh38",
            "filters": "PASS",
            "html_report_truncate_settings": {
                "columns_for_hover_expand": ["IMPACT", "ClinVar_CLNSIG", "GENE"],
                "default_max_width_px": xss_width,
            },
        }
        v_file = tmp_path / "variants.json"
        s_file = tmp_path / "summary.json"
        v_file.write_text(json.dumps(variants_data), encoding="utf-8")
        s_file.write_text(json.dumps(summary), encoding="utf-8")
        generate_html_report(str(v_file), str(s_file), str(tmp_path), cfg)
        report_html = (tmp_path / "index.html").read_text(encoding="utf-8")

        # 1. Verify onmouseover and WIDTH_PROBE are not in rendered static HTML
        assert "onmouseover" not in report_html
        assert "WIDTH_PROBE" not in report_html

        # 2. Verify in real browser that no onmouseover attribute exists
        # and WIDTH_PROBE is undefined
        probe_script = """
<script>
window.addEventListener("load", function() {
    setTimeout(function() {
        try {
            var onmouseoverEls = document.querySelectorAll("[onmouseover]");
            var hasWidthProbe = typeof globalThis.WIDTH_PROBE !== "undefined";
            var marker = document.createElement("div");
            marker.id = "width-probe-result";
            marker.textContent = JSON.stringify({
                onmouseoverCount: onmouseoverEls.length,
                hasWidthProbe: hasWidthProbe
            });
            document.body.appendChild(marker);
        } catch (e) {
            var errMarker = document.createElement("div");
            errMarker.id = "width-probe-result";
            errMarker.textContent = JSON.stringify({ error: e.message });
            document.body.appendChild(errMarker);
        }
    }, 200);
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

        server = http.server.HTTPServer(("127.0.0.1", 0), Handler)
        port = server.server_port
        t = threading.Thread(target=server.serve_forever)
        t.daemon = True
        t.start()

        try:
            cmd = [
                chrome_bin,
                "--headless",
                "--disable-gpu",
                "--no-sandbox",
                "--virtual-time-budget=2000",
                "--dump-dom",
                f"http://127.0.0.1:{port}/",
            ]
            res = subprocess.run(cmd, capture_output=True, text=True, timeout=15)
            assert res.returncode == 0, f"Headless Chrome execution failed: {res.stderr}"
            m = re.search(r'<div id="width-probe-result">([\s\S]*?)</div>', res.stdout)
            assert m is not None, f"Width probe marker not found: {res.stdout[:500]}"
            import html as html_lib

            data = json.loads(html_lib.unescape(m.group(1)))
            assert "error" not in data, f"Width probe reported error: {data.get('error')}"
            assert data["onmouseoverCount"] == 0
            assert data["hasWidthProbe"] is False
        finally:
            server.shutdown()
            server.server_close()

    def test_coerce_positive_int_overflow_and_special_values(self, tmp_path: Path):
        """Verify _coerce_positive_int handles inf, -inf, nan, and 1e400 gracefully."""
        from variantcentrifuge.generate_html_report import (
            _coerce_positive_int,
            generate_html_report,
        )

        assert _coerce_positive_int(float("inf"), default=150) == 150
        assert _coerce_positive_int(float("-inf"), default=150) == 150
        assert _coerce_positive_int(float("nan"), default=150) == 150
        assert _coerce_positive_int(1e400, default=150) == 150
        assert _coerce_positive_int("1e400", default=150) == 150

        # Verify generate_html_report does not crash with 1e400 or inf in config
        v_file = tmp_path / "variants.json"
        s_file = tmp_path / "summary.json"
        summary = {
            "num_variants": 1,
            "num_samples": 1,
            "impact_distribution": {"HIGH": 1},
            "inheritance_distribution": {"de_novo": 1},
            "top_genes": [{"gene": "BRCA1", "count": 1}],
        }
        v_file.write_text(
            json.dumps([{"CHROM": "chr1", "POS": 100, "REF": "A", "ALT": "T"}]),
            encoding="utf-8",
        )
        s_file.write_text(json.dumps(summary), encoding="utf-8")

        cfg = {
            "reference": "GRCh38",
            "html_report_truncate_settings": {
                "columns_for_hover_expand": ["CHROM"],
                "default_max_width_px": float("inf"),
            },
        }
        generate_html_report(str(v_file), str(s_file), str(tmp_path), cfg)
        out_html = (tmp_path / "index.html").read_text(encoding="utf-8")
        assert "max-width: 150px;" in out_html

    def test_symbolic_allele_igv_links_rendered_safely(self, tmp_path: Path):
        """Verify symbolic allele IGV links (<DEL>) render safely as %3CDEL%3E."""
        from variantcentrifuge.generate_html_report import (
            _encode_symbolic_allele_brackets,
            generate_html_report,
            is_safe_url,
        )

        # Direct function behavior
        assert is_safe_url("igv/S1_chr1_100_N_<DEL>_igv_report.html") is True
        assert is_safe_url("igv/S1_chr1_100_N_<DUP>_igv_report.html") is True
        assert is_safe_url("igv/S1_chr1_100_N_<INS:ME:ALU>_igv_report.html") is True
        del_path = "igv/S1_chr1_100_N_<DEL>_igv_report.html"
        expected_encoded = "igv/S1_chr1_100_N_%3CDEL%3E_igv_report.html"
        assert _encode_symbolic_allele_brackets(del_path) == expected_encoded

        # Malicious tags must be rejected
        assert is_safe_url("igv/<script>alert(1)</script>.html") is False
        assert is_safe_url("igv/<img src=x onerror=alert(1)>.html") is False
        assert is_safe_url("igv/<svg onload=alert(1)>.html") is False

        variants_data = [
            {
                "CHROM": "chr1",
                "POS": 100,
                "REF": "N",
                "ALT": "<DEL>",
                "igv_links": [
                    {
                        "sample_id": "S1_DEL",
                        "report_path": "igv/S1_chr1_100_N_<DEL>_igv_report.html",
                    },
                ],
            },
            {
                "CHROM": "chr1",
                "POS": 200,
                "REF": "A",
                "ALT": "T",
                "igv_links": [
                    {
                        "sample_id": "S2_SCRIPT",
                        "report_path": "igv/<script>alert(1)</script>.html",
                    },
                ],
            },
            {
                "CHROM": "chr1",
                "POS": 300,
                "REF": "A",
                "ALT": "T",
                "igv_links": [
                    {
                        "sample_id": "S3_IMG",
                        "report_path": "igv/<img src=x onerror=alert(1)>.html",
                    },
                ],
            },
        ]
        summary = {
            "num_variants": 3,
            "num_samples": 3,
            "impact_distribution": {"HIGH": 1},
            "inheritance_distribution": {"de_novo": 1},
            "top_genes": [{"gene": "BRCA1", "count": 1}],
        }
        v_file = tmp_path / "variants.json"
        s_file = tmp_path / "summary.json"
        v_file.write_text(json.dumps(variants_data), encoding="utf-8")
        s_file.write_text(json.dumps(summary), encoding="utf-8")

        cfg = {"reference": "GRCh38"}
        generate_html_report(str(v_file), str(s_file), str(tmp_path), cfg)
        report_html = (tmp_path / "index.html").read_text(encoding="utf-8")

        # Legitimate symbolic allele link must render as clickable with percent-encoded brackets
        assert 'href="igv/S1_chr1_100_N_%3CDEL%3E_igv_report.html"' in report_html
        assert "S1_DEL" in report_html

        # Malicious links must NOT render as clickable <a> tags
        assert 'href="igv/<script>' not in report_html
        assert 'href="igv/%3Cscript%3E' not in report_html
        assert "<script>alert(1)</script>" not in report_html
        assert 'href="igv/<img' not in report_html
        assert "S2_SCRIPT" in report_html
        assert "S3_IMG" in report_html
