"""Phase-aware compound het calls without parents (#236).

Phased GT (``0|1``/``1|0``, same PS or no PS column) or GATK PGT/PID (same PID)
prove two hets cis (pair dropped) or trans (``compound_heterozygous``,
``trans_phased``). Unphased or incomparable pairs keep the old labels; the trio
path ignores phase.
"""

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
