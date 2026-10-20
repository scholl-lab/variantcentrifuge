"""Headless-Chrome layout regression for narrow viewports (#185)."""

import html
import http.server
import json
import re
import shutil
import subprocess
import threading
from pathlib import Path

import pytest

from variantcentrifuge.generate_html_report import generate_html_report

PROBE = """
<script>
window.__pageErrors = 0;
window.addEventListener("error", function() { window.__pageErrors++; });
window.addEventListener("load", function() {
    setTimeout(function() {
        function measure() {
            var vw = window.innerWidth;
            var right = 0;
            document.querySelectorAll(
                "#table-toolbar button, #table-toolbar input, #table-toolbar select, " +
                "#table-toolbar label, #filter-row .filter-control"
            ).forEach(function(el) {
                var r = el.getBoundingClientRect();
                if (r.width > 0) right = Math.max(right, Math.ceil(r.right));
            });
            var scroller = document.querySelector(".dt-scroll-body");
            return {
                viewport: vw,
                document: document.documentElement.scrollWidth,
                controlsRight: right,
                tableScrollWidth: scroller ? scroller.scrollWidth : -1,
                tableClientWidth: scroller ? scroller.clientWidth : -1
            };
        }
        var out = {};
        try {
            out.collapsed = measure();
            document.getElementById("filter-toggle-btn").click();
            out.filterControls = document.querySelectorAll("#filter-row .filter-control").length;
            out.expanded = measure();
            out.filterRowShown = document.getElementById("filter-row").offsetHeight > 0;
            out.toolbarHeight = document.getElementById("table-toolbar").offsetHeight;
            out.rows = $("#variants_table").DataTable().rows().count();
            out.external = performance.getEntriesByType("resource").filter(function(e) {
                return e.name.indexOf(location.origin) !== 0;
            }).length;
        } catch (e) {
            out.error = e.message;
        }
        out.pageErrors = window.__pageErrors;
        // The report runs in a fixed-size iframe; report back to the dumped top document.
        var marker = window.parent.document.createElement("div");
        marker.id = "layout-probe-result";
        marker.textContent = JSON.stringify(out);
        window.parent.document.body.appendChild(marker);
    }, 500);
});
</script>
"""


def _render_report(tmp_path: Path) -> str:
    # Wide synthetic table: many columns, several numeric ones so the filter panel fills.
    variants = [
        {
            "CHROM": f"chr{i % 22 + 1}",
            "POS": 1000 + i,
            "REF": "A",
            "ALT": "G",
            "GENE": f"GENE{i}",
            "IMPACT": ["HIGH", "MODERATE", "LOW"][i % 3],
            "EFFECT": "missense_variant",
            "HGVS_C": f"c.{100 + i}A>G",
            "HGVS_P": f"p.Lys{30 + i}Arg",
            "CADD_PHRED": 10.0 + i,
            "REVEL": round(0.05 * i, 3),
            "gnomAD_AF": 1e-4 * (i + 1),
            "AC": i + 1,
            "DP": 30 + i,
            "GQ": 60 + i,
            "nomen_long_annotation_column": "x" * 40,
            # Long categorical values must not widen the filter panel (#185 review).
            "ClinVar_CLNSIG": (
                "Conflicting_classifications_of_pathogenicity|Pathogenic|"
                "Likely_pathogenic|Uncertain_significance"
                if i % 2
                else "Benign"
            ),
        }
        for i in range(12)
    ]
    summary = {
        "num_variants": 12,
        "num_genes": 12,
        "num_samples": 1,
        "impact_distribution": {"HIGH": 4, "MODERATE": 4, "LOW": 4},
    }
    v_file = tmp_path / "variants.json"
    s_file = tmp_path / "summary.json"
    v_file.write_text(json.dumps(variants), encoding="utf-8")
    s_file.write_text(json.dumps(summary), encoding="utf-8")
    # Realistic footer metadata: a long filter expression and VCF name must not widen the page.
    cfg = {
        "filters": "((ANN[0].IMPACT has 'HIGH') | (ANN[0].IMPACT has 'MODERATE'))",
        "vcf_file": "/data/synthetic_cohort_panel.GRCh38.annotated.normalized.vcf.gz",
        "reference": "GRCh38.p14",
    }
    generate_html_report(str(v_file), str(s_file), str(tmp_path), cfg)
    return (tmp_path / "index.html").read_text(encoding="utf-8")


def _probe_layout(chrome_bin: str, page: str, width: int, height: int) -> dict:
    # Headless Chrome clamps --window-size to >= 500px, so the viewport is an iframe.
    frame = (
        "<!doctype html><body style='margin:0'><iframe src='/report' "
        f"style='border:0;width:{width}px;height:{height}px'></iframe></body>"
    )

    class Handler(http.server.BaseHTTPRequestHandler):
        def do_GET(self):
            self.send_response(200)
            self.send_header("Content-Type", "text/html; charset=utf-8")
            self.end_headers()
            self.wfile.write((page if self.path == "/report" else frame).encode("utf-8"))

        def log_message(self, *args):
            pass

    server = http.server.ThreadingHTTPServer(("127.0.0.1", 0), Handler)
    threading.Thread(target=server.serve_forever, daemon=True).start()
    cmd = [
        chrome_bin,
        "--headless",
        "--disable-gpu",
        "--no-sandbox",
        "--no-first-run",
        "--hide-scrollbars",
        f"--window-size={max(width, 500) + 50},{height + 50}",
        "--virtual-time-budget=6000",
        "--dump-dom",
        f"http://127.0.0.1:{server.server_port}/",
    ]
    try:
        try:
            res = subprocess.run(cmd, capture_output=True, text=True, timeout=45)
        except subprocess.TimeoutExpired:  # headless Chrome occasionally hangs on loaded CI
            res = subprocess.run(cmd, capture_output=True, text=True, timeout=45)
    finally:
        server.shutdown()
        server.server_close()
    assert res.returncode == 0, f"Chrome failed: {res.stderr}"
    m = re.search(r'<div id="layout-probe-result">([\s\S]*?)</div>', res.stdout)
    assert m is not None, f"Probe marker not found: {res.stdout[-500:]}"
    data = json.loads(html.unescape(m.group(1)))
    assert "error" not in data, data
    return data


@pytest.mark.unit
@pytest.mark.parametrize("width,height", [(390, 844), (1440, 900)])
def test_toolbar_and_filters_fit_viewport_with_scrollable_table(tmp_path, width, height):
    chrome_bin = (
        shutil.which("google-chrome")
        or shutil.which("chromium")
        or shutil.which("chromium-browser")
    )
    if not chrome_bin:
        pytest.skip("Chrome/Chromium binary not available")

    page = _render_report(tmp_path).replace("</body>", PROBE + "</body>")
    data = _probe_layout(chrome_bin, page, width, height)

    assert data["collapsed"]["viewport"] == width
    assert data["rows"] == 12
    assert data["filterControls"] > 0
    assert data["filterRowShown"] is True
    assert data["pageErrors"] == 0
    assert data["external"] == 0
    for state in ("collapsed", "expanded"):
        m = data[state]
        assert m["document"] <= m["viewport"], (state, m)
        assert m["controlsRight"] <= m["viewport"], (state, m)
    if width >= 1024:
        assert data["toolbarHeight"] < 40, data  # desktop toolbar stays a single row
    else:
        # The wide variants table keeps its own horizontal scroll container.
        m = data["collapsed"]
        assert m["tableScrollWidth"] > m["tableClientWidth"] > 0, m
