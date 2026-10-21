"""Phase-aware compound het calls without parents (#236).

Phased GT (``0|1``/``1|0``, same PS or no PS column) or GATK PGT/PID (same PID)
prove two hets cis (pair dropped) or trans (``compound_heterozygous``,
``trans_phased``). Unphased or incomparable pairs keep the old labels; the trio
path ignores phase.
"""

import logging
from unittest.mock import patch

import pandas as pd

from variantcentrifuge.inheritance.analyzer import analyze_inheritance
from variantcentrifuge.inheritance.comp_het_vectorized import (
    analyze_gene_for_compound_het_vectorized,
)

K1, K2, K3 = "1:100:A>T", "1:200:C>G", "1:300:G>A"
NO_PED = "compound_heterozygous_possible_no_pedigree"


def _gene(child_gts, **extra):
    rows = []
    for i, gt in enumerate(child_gts):
        row = {
            "CHROM": "1",
            "POS": str(100 * (i + 1)),
            "REF": "ACG"[i],
            "ALT": "TGA"[i],
            "GENE": "G1",
            "child": gt,
        }
        for col, vals in extra.items():
            row[col] = vals[i]
        rows.append(row)
    return pd.DataFrame(rows)


def _run(df, samples=("child",), ped=None):
    return analyze_gene_for_compound_het_vectorized(df, ped or {}, list(samples))


def test_cis_pair_is_dropped():
    assert _run(_gene(["0|1", "0|1"])) == {}
    assert _run(_gene(["1|0", "1|0"])) == {}


def test_trans_pair_is_compound_het():
    res = _run(_gene(["0|1", "1|0"]))
    for key, partner in ((K1, K2), (K2, K1)):
        info = res[key]["child"]
        assert info["comp_het_type"] == "compound_heterozygous"
        assert info["inheritance_type"] == "trans_phased"
        assert info["partner_variants"] == [partner]
        assert info["partner_types"] == {partner: "trans_phased"}
        assert info["partner_comp_het_types"] == {partner: "compound_heterozygous"}
        assert info["all_comp_het_types"] == ["compound_heterozygous"]


def test_three_hets_mixed_cis_and_trans():
    res = _run(_gene(["0|1", "0|1", "1|0"]))
    assert res[K1]["child"]["partner_variants"] == [K3]
    assert res[K2]["child"]["partner_variants"] == [K3]
    assert res[K3]["child"]["partner_variants"] == [K1, K2]
    assert {res[k]["child"]["comp_het_type"] for k in (K1, K2, K3)} == {"compound_heterozygous"}


def test_phased_and_unphased_mixed():
    res = _run(_gene(["0|1", "1|0", "0/1"]))
    info = res[K1]["child"]
    assert info["comp_het_type"] == "compound_heterozygous"
    assert info["inheritance_type"] == "trans_phased"
    assert info["partner_variants"] == [K2, K3]
    assert info["partner_types"] == {K2: "trans_phased", K3: "unknown"}
    assert info["partner_comp_het_types"] == {K2: "compound_heterozygous", K3: NO_PED}
    assert info["all_comp_het_types"] == ["compound_heterozygous", NO_PED]
    info3 = res[K3]["child"]
    assert (info3["comp_het_type"], info3["inheritance_type"]) == (NO_PED, "unknown")
    assert info3["partner_variants"] == [K1, K2]


def test_phase_set_must_match():
    # different PS -> phase not comparable -> old label
    res = _run(_gene(["0|1", "0|1"], **{"GEN[0].PS": ["100", "200"]}))
    assert res[K1]["child"]["comp_het_type"] == NO_PED
    # same PS (sanitized column name) -> cis -> dropped
    assert _run(_gene(["0|1", "0|1"], **{"GEN_0__PS": ["100", "100"]})) == {}
    # PS column present but missing for one record -> unknown
    res = _run(_gene(["0|1", "1|0"], **{"GEN[0].PS": ["100", "."]}))
    assert res[K1]["child"]["comp_het_type"] == NO_PED


def test_ps_column_of_other_sample_is_ignored():
    df = _gene(["0|1", "1|0"], **{"GEN[1].PS": ["100", "200"], "other": ["0/0", "0/0"]})
    res = _run(df, samples=("child", "other"))
    assert res[K1]["child"]["comp_het_type"] == "compound_heterozygous"


def test_pgt_pid():
    def gene(pgt, pid):
        return _gene(["0/1", "0/1"], **{"GEN[0].PGT": pgt, "GEN[0].PID": pid})

    res = _run(gene(["0|1", "1|0"], ["100_A_T", "100_A_T"]))
    assert res[K1]["child"]["comp_het_type"] == "compound_heterozygous"
    assert res[K1]["child"]["inheritance_type"] == "trans_phased"
    assert _run(gene(["0|1", "0|1"], ["100_A_T", "100_A_T"])) == {}
    res = _run(gene(["0|1", "0|1"], ["100_A_T", "300_G_A"]))
    assert res[K1]["child"]["comp_het_type"] == NO_PED
    res = _run(gene(["0|1", "."], ["100_A_T", "."]))
    assert res[K1]["child"]["comp_het_type"] == NO_PED


def test_unphased_labels_unchanged():
    res = _run(_gene(["0/1", "1/0"]))
    assert res[K1]["child"] == {
        "sample_id": "child",
        "gene": "G1",
        "partner_variants": [K2],
        "is_compound_het": True,
        "comp_het_type": NO_PED,
        "inheritance_type": "unknown",
        "partner_types": {K2: "unknown"},
        "partner_comp_het_types": {K2: NO_PED},
        "all_comp_het_types": [NO_PED],
    }


def _ped(father="0", mother="0"):
    return {
        "child": {
            "sample_id": "child",
            "father_id": father,
            "mother_id": mother,
            "affected_status": "2",
        },
        "father": {
            "sample_id": "father",
            "father_id": "0",
            "mother_id": "0",
            "affected_status": "1",
        },
        "mother": {
            "sample_id": "mother",
            "father_id": "0",
            "mother_id": "0",
            "affected_status": "1",
        },
    }


def test_single_parent_cis_dropped_trans_called():
    samples = ("child", "father")
    base = {"father": ["0|1", "0|0"]}
    cis = _run(_gene(["0|1", "0|1"], **base), samples, _ped(father="father"))
    assert cis == {}
    res = _run(_gene(["0|1", "1|0"], **base), samples, _ped(father="father"))
    assert res[K1]["child"]["comp_het_type"] == "compound_heterozygous"
    unph = _run(_gene(["0/1", "0/1"], father=["0/1", "0/0"]), samples, _ped(father="father"))
    assert unph[K1]["child"]["comp_het_type"] == "compound_heterozygous_possible_missing_parents"


def test_trio_ignores_phase():
    # parents say trans; child's phased GT says cis -> trio result stands
    df = _gene(["0|1", "0|1"], father=["0|1", "0|0"], mother=["0|0", "0|1"])
    res = _run(df, ("child", "father", "mother"), _ped("father", "mother"))
    assert res[K1]["child"]["comp_het_type"] == "compound_heterozygous"
    assert res[K1]["child"]["inheritance_type"] == "trans"


def test_end_to_end_patterns():
    df = _gene(["0|1", "1|0", "0|0"])
    out = analyze_inheritance(df, {}, ["child"])
    assert out["Inheritance_Pattern"].tolist()[:2] == ["compound_heterozygous"] * 2
    out = analyze_inheritance(_gene(["0|1", "0|1"]), {}, ["child"])
    assert NO_PED not in out["Inheritance_Pattern"].tolist()


# --- review r1: phase blocks declared in the header but not extracted -------------
def _hdr(*format_ids):
    return "\n".join(
        f'##FORMAT=<ID={i},Number=1,Type=String,Description="x">' for i in ("GT", *format_ids)
    )


def _usable(header, fields):
    from variantcentrifuge.vcf_header_parser import comp_het_phase_usable

    cfg = {"vcf_file": "in.vcf.gz", "fields_to_extract": fields}
    with patch("variantcentrifuge.vcf_header_parser._read_phase_header", return_value=header):
        return comp_het_phase_usable(cfg), cfg


def test_header_without_ps_or_pid_uses_phase():
    assert _usable(_hdr(), "CHROM POS GEN[*].GT")[0] is True


def test_header_ps_not_extracted_turns_phase_off(caplog):
    with caplog.at_level(logging.INFO, logger="variantcentrifuge"):
        ok, cfg = _usable(_hdr("PS"), "CHROM POS GEN[*].GT")
    assert ok is False
    msgs = [r.getMessage() for r in caplog.records if "phase" in r.getMessage()]
    assert len(msgs) == 1 and "GEN[*].PS" in msgs[0]
    # cached: second call neither re-reads the header nor logs again
    with patch("variantcentrifuge.vcf_header_parser._read_phase_header") as hdr:
        from variantcentrifuge.vcf_header_parser import comp_het_phase_usable

        assert comp_het_phase_usable(cfg) is False
        hdr.assert_not_called()


def test_unreadable_header_fails_closed(caplog):
    import subprocess

    from variantcentrifuge.vcf_header_parser import comp_het_phase_usable

    err = subprocess.CalledProcessError(255, ["bcftools"], stderr="cannot open")
    cfg = {"vcf_file": "missing.vcf.gz", "fields_to_extract": "CHROM GEN[*].GT"}
    with (
        caplog.at_level(logging.INFO, logger="variantcentrifuge"),
        patch("variantcentrifuge.vcf_header_parser._read_phase_header", side_effect=err),
    ):
        assert comp_het_phase_usable(cfg) is False
    assert not [r for r in caplog.records if r.levelno >= logging.ERROR]
    assert any(r.levelno == logging.WARNING and "phase" in r.getMessage() for r in caplog.records)
    # PS extracted: a missing header no longer matters for PS blocks
    cfg = {"vcf_file": "missing.vcf.gz", "fields_to_extract": "GEN[*].GT GEN[*].PS"}
    with patch("variantcentrifuge.vcf_header_parser._read_phase_header", side_effect=err):
        assert comp_het_phase_usable(cfg) is True


def test_compound_het_description_unchanged():
    from variantcentrifuge.inheritance.prioritizer import get_pattern_description

    assert get_pattern_description("compound_heterozygous") == (
        "Two different mutations in the same gene, one from each parent"
    )


def test_header_pid_not_extracted_turns_phase_off():
    assert _usable(_hdr("PGT", "PID"), "CHROM GEN[*].GT GEN[*].PS")[0] is False


def test_header_ps_extracted_uses_phase():
    assert _usable(_hdr("PS"), "CHROM POS GEN[*].GT GEN[*].PS")[0] is True
    assert _usable(_hdr("PS", "PGT", "PID"), "GEN[*].GT GEN[*].PS GEN[*].PGT GEN[*].PID")[0]


def test_use_phase_false_keeps_old_labels():
    res = analyze_gene_for_compound_het_vectorized(
        _gene(["0|1", "0|1"]), {}, ["child"], use_phase=False
    )
    assert res[K1]["child"]["comp_het_type"] == NO_PED
    out = analyze_inheritance(_gene(["0|1", "1|0"]), {}, ["child"], use_phase=False)
    assert out["Inheritance_Pattern"].tolist() == [NO_PED] * 2


def test_parallel_analyzer_passes_use_phase():
    from variantcentrifuge.inheritance.parallel_analyzer import analyze_inheritance_parallel

    df = pd.concat([_gene(["0|1", "0|1"])] * 1, ignore_index=True)
    out = analyze_inheritance_parallel(
        df, {}, ["child"], n_workers=2, min_variants_for_parallel=1, use_phase=False
    )
    assert out["Inheritance_Pattern"].tolist() == [NO_PED] * 2
    out = analyze_inheritance_parallel(df, {}, ["child"], n_workers=2, min_variants_for_parallel=1)
    assert NO_PED not in out["Inheritance_Pattern"].tolist()


def test_stage_reads_header_and_passes_use_phase():
    from tests.mocks.fixtures import create_test_context
    from variantcentrifuge.stages.analysis_stages import InheritanceAnalysisStage

    ctx = create_test_context()
    ctx.current_dataframe = _gene(["0|1", "0|1"])
    ctx.vcf_samples = ["child"]
    ctx.pedigree_data = {}
    ctx.config.update(
        inheritance_mode="simple", vcf_file="in.vcf.gz", fields_to_extract="CHROM GEN[*].GT"
    )
    for s in ("dataframe_loading", "sample_config_loading", "pedigree_loading"):
        ctx.mark_complete(s)
    with patch("variantcentrifuge.vcf_header_parser._read_phase_header", return_value=_hdr("PS")):
        out = InheritanceAnalysisStage()._process(ctx)
    assert out.current_dataframe["Inheritance_Pattern"].tolist() == [NO_PED] * 2


def test_chunk_path_passes_use_phase():
    from tests.mocks.fixtures import create_test_context
    from variantcentrifuge.stages.chunk_helpers import apply_chunk_inheritance_analysis

    ctx = create_test_context()
    ctx.config["inheritance_analysis_config"] = {
        "vcf_samples": ["child"],
        "pedigree_data": {},
        "use_phase": False,
    }
    out = apply_chunk_inheritance_analysis(_gene(["0|1", "0|1"]), ctx)
    assert out["Inheritance_Pattern"].tolist() == [NO_PED] * 2


# --- whole-branch review: header key order, PID without PGT, effective fields, no VCF ---
REORDERED = (
    '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'
    '##FORMAT=<ID=PS,Type=Integer,Number=1,Description="Phase set">\n'
    '##FORMAT=<Type=String,ID=PID,Number=1,Description="Physical phasing ID, a, b">\n'
    '##INFO=<Description="Allele frequency",Number=A,Type=Float,ID=AF>\n'
)


def test_reordered_header_keys_are_parsed():
    from variantcentrifuge.vcf_header_parser import parse_vcf_header

    with patch("variantcentrifuge.vcf_header_parser._read_header", return_value=REORDERED):
        got = {(f.source, f.id, f.number, f.type, f.description) for f in parse_vcf_header("x")}
    assert ("FORMAT", "PS", "1", "Integer", "Phase set") in got
    assert ("FORMAT", "PID", "1", "String", "Physical phasing ID, a, b") in got
    assert ("INFO", "AF", "A", "Float", "Allele frequency") in got


def test_reordered_ps_header_turns_phase_off():
    assert _usable(REORDERED, "CHROM GEN[*].GT")[0] is False
    assert _usable(REORDERED, "CHROM GEN[*].GT GEN[*].PS GEN[*].PID")[0] is True


def test_pid_without_pgt_never_uses_global_block():
    # PID extracted, PGT not: phased GT without PS must stay unknown (old label)
    res = _run(_gene(["0|1", "0|1"], **{"GEN[0].PID": ["100_A_T", "300_G_A"]}))
    assert res[K1]["child"]["comp_het_type"] == NO_PED
    # PGT present but missing for these records, no PS: unknown as well
    df = _gene(["0|1", "1|0"], **{"GEN[0].PID": [".", "."], "GEN[0].PGT": [".", "."]})
    assert _run(df)[K1]["child"]["comp_het_type"] == NO_PED
    # GT with an extracted PS still proves phase even when PID is present
    df = _gene(["0|1", "1|0"], **{"GEN[0].PID": [".", "."], "GEN[0].PS": ["7", "7"]})
    assert _run(df)[K1]["child"]["comp_het_type"] == "compound_heterozygous"


def test_appended_extra_sample_fields_count_as_extracted():
    from variantcentrifuge.vcf_header_parser import comp_het_phase_usable

    cfg = {
        "vcf_file": "in.vcf.gz",
        "fields_to_extract": "CHROM GEN[*].GT",
        "append_extra_sample_fields": True,
        "extra_sample_fields": ["PS"],
    }
    with patch("variantcentrifuge.vcf_header_parser._read_phase_header", return_value=_hdr("PS")):
        assert comp_het_phase_usable(cfg) is True


def test_no_vcf_returns_false_without_bcftools_or_warning(caplog):
    from variantcentrifuge.vcf_header_parser import comp_het_phase_usable

    with (
        caplog.at_level(logging.INFO, logger="variantcentrifuge"),
        patch("variantcentrifuge.vcf_header_parser._read_phase_header") as rd,
    ):
        assert comp_het_phase_usable({"fields_to_extract": "GEN[*].GT GEN[*].PS"}) is False
    rd.assert_not_called()
    assert not [r for r in caplog.records if r.levelno >= logging.WARNING]


def test_analyze_variants_passes_use_phase():
    import variantcentrifuge.analyze_variants as av

    tsv = "CHROM\tPOS\tREF\tALT\tGENE\tGT\n1\t100\tA\tT\tG1\tchild(0|1)\n"
    cfg = {
        "calculate_inheritance": True,
        "pedigree_data": {},
        "sample_list": "child",
        "no_stats": True,
    }
    with (
        patch(
            "variantcentrifuge.inheritance.analyze_inheritance", side_effect=lambda df, *a, **k: df
        ) as ai,
        patch.object(av, "comp_het_phase_usable", return_value=False) as usable,
    ):
        list(av.analyze_variants(iter([tsv]), cfg))
    usable.assert_called_once_with(cfg)
    assert ai.call_args.kwargs["use_phase"] is False


# --- B3b codex re-review: attribute boundaries, exact field tokens, quote pair ---------
DOTTED_AFTER = (
    '##FORMAT=<ID=PS,Number=1,Type=Integer,Description="Phase set, cutoff=3 > 1",'
    "Original.ID=phase_set>\n"
    '##FORMAT=<ID=PID,Number=1,Type=String,Description="Phasing ID",Original.ID=pid_x>\n'
    '##INFO=<ID=AF,Number=A,Type=Float,Description="AF",Original.ID=af_raw>\n'
)
DOTTED_BEFORE = (
    '##FORMAT=<Original.ID=phase_set,ID=PS,Number=1,Type=Integer,Description="Phase set">\n'
    '##FORMAT=<Original.ID=pid_x,ID=PID,Number=1,Type=String,Description="Phasing ID">\n'
    '##INFO=<Original.ID=af_raw,ID=AF,Number=A,Type=Float,Description="AF">\n'
)


def _parsed(text):
    from variantcentrifuge.vcf_header_parser import parse_vcf_header

    with patch("variantcentrifuge.vcf_header_parser._read_header", return_value=text):
        return {(f.source, f.id): f for f in parse_vcf_header("x")}


def test_dotted_attributes_do_not_overwrite_id():
    for text in (DOTTED_AFTER, DOTTED_BEFORE):
        got = _parsed(text)
        assert {("FORMAT", "PS"), ("FORMAT", "PID"), ("INFO", "AF")} == set(got)
        assert _usable(text, "CHROM GEN[*].GT")[0] is False
    assert _parsed(DOTTED_AFTER)[("FORMAT", "PS")].description == "Phase set, cutoff=3 > 1"


def test_description_keeps_escaped_quotes():
    text = (
        '##INFO=<ID=X,Number=1,Type=String,Description="a, b=c > d said \\"PASS\\"">\n'
        '##INFO=<ID=Y,Number=1,Type=String,Description="\\"q\\" first">\n'
    )
    got = _parsed(text)
    assert got[("INFO", "X")].description == 'a, b=c > d said \\"PASS\\"'
    assert got[("INFO", "Y")].description == '\\"q\\" first'


def test_appended_ps_aux_does_not_count_as_ps():
    from variantcentrifuge.vcf_header_parser import comp_het_phase_usable

    cfg = {
        "vcf_file": "in.vcf.gz",
        "fields_to_extract": "CHROM GEN[*].GT",
        "append_extra_sample_fields": True,
        "extra_sample_fields": ["GEN[*].PS.AUX"],
    }
    with patch("variantcentrifuge.vcf_header_parser._read_phase_header", return_value=_hdr("PS")):
        use_phase = comp_het_phase_usable(cfg)
    assert use_phase is False
    df = _gene(["0|1", "0|1"], **{"GEN[0].PS.AUX": ["100", "200"]})
    res = analyze_gene_for_compound_het_vectorized(df, {}, ["child"], use_phase=use_phase)
    assert res[K1]["child"]["comp_het_type"] == NO_PED
