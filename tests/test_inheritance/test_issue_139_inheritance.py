"""
Regression tests for Issue #139:
1. Reject de_novo_candidate if a non-missing called parent carries ALT.
2. Restrict autosomal_dominant and autosomal_recessive patterns to autosomes only.
3. Make from_neither and from_both reciprocal in compound-het partner detection.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from variantcentrifuge.inheritance.comp_het_vectorized import (
    find_potential_partners_vectorized,
)
from variantcentrifuge.inheritance.deducer import (
    check_dominant_pattern,
    check_recessive_pattern,
    deduce_patterns_for_sample,
)
from variantcentrifuge.inheritance.vectorized_deducer import (
    vectorized_deduce_patterns,
)


@pytest.fixture
def trio_pedigree() -> dict[str, dict[str, str]]:
    return {
        "proband": {
            "sample_id": "proband",
            "father_id": "father",
            "mother_id": "mother",
            "sex": "1",  # male
            "affected_status": "2",  # affected
        },
        "father": {
            "sample_id": "father",
            "father_id": "0",
            "mother_id": "0",
            "sex": "1",
            "affected_status": "1",  # unaffected
        },
        "mother": {
            "sample_id": "mother",
            "father_id": "0",
            "mother_id": "0",
            "sex": "2",
            "affected_status": "1",  # unaffected
        },
    }


@pytest.mark.unit
def test_de_novo_candidate_rejected_when_parent_carries_alt(trio_pedigree):
    """When a parent carries ALT, the variant is inherited, not a de novo candidate."""
    # Case 1: Father is 0/1 (carrier), Mother is ./. (missing)
    row_inherited = {
        "CHROM": "chr1",
        "POS": 100,
        "proband": "0/1",
        "father": "0/1",
        "mother": "./.",
    }
    sample_list = ["proband", "father", "mother"]

    patterns_deducer = deduce_patterns_for_sample(
        "proband", row_inherited, trio_pedigree, sample_list
    )
    assert "de_novo_candidate" not in patterns_deducer
    assert "de_novo" not in patterns_deducer

    df_inherited = pd.DataFrame([row_inherited])
    patterns_vec = vectorized_deduce_patterns(df_inherited, trio_pedigree, sample_list)
    assert "de_novo_candidate" not in patterns_vec[0]
    assert "de_novo" not in patterns_vec[0]

    # Case 2: Father is 0/0 (reference), Mother is ./. (missing) -> valid de_novo_candidate
    row_candidate = {
        "CHROM": "chr1",
        "POS": 200,
        "proband": "0/1",
        "father": "0/0",
        "mother": "./.",
    }
    patterns_deducer_cand = deduce_patterns_for_sample(
        "proband", row_candidate, trio_pedigree, sample_list
    )
    assert "de_novo_candidate" in patterns_deducer_cand

    df_candidate = pd.DataFrame([row_candidate])
    patterns_vec_cand = vectorized_deduce_patterns(df_candidate, trio_pedigree, sample_list)
    assert "de_novo_candidate" in patterns_vec_cand[0]


@pytest.mark.unit
def test_autosomal_patterns_restricted_to_autosomes(trio_pedigree):
    """X, Y, M, MT variants must never receive autosomal_* patterns."""
    sample_list = ["proband", "father", "mother"]
    non_autosomal_chroms = ["X", "chrX", "Y", "chrY", "MT", "chrM", "23"]

    for chrom in non_autosomal_chroms:
        row = {
            "CHROM": chrom,
            "POS": 100,
            "proband": "0/1",
            "father": "0/1",
            "mother": "0/0",
        }
        # Direct helper check
        assert check_dominant_pattern("proband", row, trio_pedigree, sample_list) is None
        assert check_recessive_pattern("proband", row, trio_pedigree, sample_list) is None

        # deducer check
        pats = deduce_patterns_for_sample("proband", row, trio_pedigree, sample_list)
        for p in pats:
            assert not p.startswith("autosomal_"), f"Unexpected pattern {p} for chrom {chrom}"

        # vectorized deducer check
        df = pd.DataFrame([row])
        vec_pats = vectorized_deduce_patterns(df, trio_pedigree, sample_list)
        for p in vec_pats[0]:
            assert not p.startswith("autosomal_"), (
                f"Unexpected vectorized pattern {p} for chrom {chrom}"
            )


@pytest.mark.unit
def test_reciprocal_compound_het_partners():
    """Verify B in partners[A] iff A in partners[B] for all variants, with no self-links."""
    # 4 heterozygous variants in proband:
    # 0: father 0/1, mother 0/0 (from father only)
    # 1: father 0/0, mother 0/1 (from mother only)
    # 2: father 0/1, mother 0/1 (from both)
    # 3: father 0/0, mother 0/0 (from neither / de novo / missing)
    father_gts = np.array([1, 0, 1, 0])
    mother_gts = np.array([0, 1, 1, 0])

    het_indices = np.array([0, 1, 2, 3])
    partners = find_potential_partners_vectorized(het_indices, father_gts, mother_gts)

    # 1. No self-links
    for idx in range(4):
        assert idx not in partners.get(idx, []), f"Variant {idx} is partnered with itself"

    # 2. Perfect reciprocity: j in partners[i] iff i in partners[j]
    for i in range(4):
        for j in range(4):
            if i == j:
                continue
            i_has_j = j in partners.get(i, [])
            j_has_i = i in partners.get(j, [])
            assert i_has_j == j_has_i, (
                f"Asymmetry between {i} and {j}: {i} has {j} = {i_has_j}, "
                f"but {j} has {i} = {j_has_i}"
            )

    # 3. Specific expected reciprocal partnerships:
    # Father-only (0) and mother-only (1) must pair
    assert 1 in partners[0] and 0 in partners[1]
    # Father-only (0) and from_both (2) must pair
    assert 2 in partners[0] and 0 in partners[2]
    # Father-only (0) and from_neither (3) must pair
    assert 3 in partners[0] and 0 in partners[3]
    # Mother-only (1) and from_both (2) must pair
    assert 2 in partners[1] and 1 in partners[2]
    # Mother-only (1) and from_neither (3) must pair
    assert 3 in partners[1] and 1 in partners[3]
    # from_both (2) and from_neither (3) both lack clear origin so they do not pair
    assert 3 not in partners.get(2, []) and 2 not in partners.get(3, [])


@pytest.mark.unit
def test_reciprocal_partner_permutations_24():
    """Verify reciprocity holds across all 24 permutations of the 4 origin classes."""
    import itertools

    base_cases = [
        (1, 0),  # 0: father-only
        (0, 1),  # 1: mother-only
        (1, 1),  # 2: from-both
        (0, 0),  # 3: from-neither
    ]

    for perm in itertools.permutations(range(4)):
        father_gts = np.array([base_cases[i][0] for i in perm])
        mother_gts = np.array([base_cases[i][1] for i in perm])
        het_indices = np.arange(4)

        partners = find_potential_partners_vectorized(het_indices, father_gts, mother_gts)

        for i in range(4):
            assert i not in partners.get(i, []), f"Self-partner in perm {perm}"
            for j in range(4):
                if i == j:
                    continue
                orig_i = perm[i]
                orig_j = perm[j]
                i_has_j = j in partners.get(i, [])
                j_has_i = i in partners.get(j, [])
                assert i_has_j == j_has_i, (
                    f"Permutation {perm}: asymmetry between pos {i} (orig {orig_i}) "
                    f"and pos {j} (orig {orig_j})"
                )


@pytest.mark.unit
def test_compound_het_ambiguous_vs_confirmed_trans(trio_pedigree):
    """Verify confirmed trans vs ambiguous compound het classification."""
    from variantcentrifuge.inheritance.comp_het_vectorized import (
        analyze_gene_for_compound_het_vectorized,
    )

    sample_list = ["proband", "father", "mother"]

    # 1. Confirmed trans: var1 father-only, var2 mother-only
    df_trans = pd.DataFrame(
        [
            {
                "CHROM": "chr1",
                "POS": "100",
                "REF": "A",
                "ALT": "T",
                "GENE": "TEST_GENE",
                "proband": "0/1",
                "father": "0/1",
                "mother": "0/0",
            },
            {
                "CHROM": "chr1",
                "POS": "200",
                "REF": "C",
                "ALT": "G",
                "GENE": "TEST_GENE",
                "proband": "0/1",
                "father": "0/0",
                "mother": "0/1",
            },
        ]
    )
    res_trans = analyze_gene_for_compound_het_vectorized(df_trans, trio_pedigree, sample_list)
    assert len(res_trans) == 2
    for _var_key, sample_dict in res_trans.items():
        assert sample_dict["proband"]["comp_het_type"] == "compound_heterozygous"
        assert sample_dict["proband"]["inheritance_type"] == "trans"

    # 2. Ambiguous: var1 father-only, var2 from both
    df_ambig = pd.DataFrame(
        [
            {
                "CHROM": "chr1",
                "POS": "100",
                "REF": "A",
                "ALT": "T",
                "GENE": "TEST_GENE",
                "proband": "0/1",
                "father": "0/1",
                "mother": "0/0",
            },
            {
                "CHROM": "chr1",
                "POS": "200",
                "REF": "C",
                "ALT": "G",
                "GENE": "TEST_GENE",
                "proband": "0/1",
                "father": "0/1",
                "mother": "0/1",
            },
        ]
    )
    res_ambig = analyze_gene_for_compound_het_vectorized(df_ambig, trio_pedigree, sample_list)
    assert len(res_ambig) == 2
    for _var_key, sample_dict in res_ambig.items():
        assert sample_dict["proband"]["comp_het_type"] == "compound_heterozygous_possible"
        assert sample_dict["proband"]["inheritance_type"] == "ambiguous"

    # 3. Missing parent genotypes: var1 father-only, var2 mother missing
    df_missing = pd.DataFrame(
        [
            {
                "CHROM": "chr1",
                "POS": "100",
                "REF": "A",
                "ALT": "T",
                "GENE": "TEST_GENE",
                "proband": "0/1",
                "father": "0/1",
                "mother": "0/0",
            },
            {
                "CHROM": "chr1",
                "POS": "200",
                "REF": "C",
                "ALT": "G",
                "GENE": "TEST_GENE",
                "proband": "0/1",
                "father": "0/0",
                "mother": "./.",
            },
        ]
    )
    res_missing = analyze_gene_for_compound_het_vectorized(df_missing, trio_pedigree, sample_list)
    assert len(res_missing) == 2
    for _var_key, sample_dict in res_missing.items():
        assert (
            sample_dict["proband"]["comp_het_type"]
            == "compound_heterozygous_possible_missing_parent_genotypes"
        )
        assert sample_dict["proband"]["inheritance_type"] == "unknown"


@pytest.mark.unit
def test_de_novo_comprehensive_parental_combinations(trio_pedigree):
    """Verify de_novo and de_novo_candidate under all parental combinations."""
    sample_list = ["proband", "father", "mother"]

    # 1. 0/0 + 0/0 -> confirmed de_novo
    row_denovo = {"CHROM": "chr1", "POS": 100, "proband": "0/1", "father": "0/0", "mother": "0/0"}
    pats = deduce_patterns_for_sample("proband", row_denovo, trio_pedigree, sample_list)
    assert "de_novo" in pats
    assert "de_novo_candidate" not in pats

    # 2. 0/0 + ./. -> de_novo_candidate
    row_fref_mmiss = {
        "CHROM": "chr1",
        "POS": 100,
        "proband": "0/1",
        "father": "0/0",
        "mother": "./.",
    }
    pats = deduce_patterns_for_sample("proband", row_fref_mmiss, trio_pedigree, sample_list)
    assert "de_novo_candidate" in pats
    assert "de_novo" not in pats

    # 3. ./. + 0/0 -> de_novo_candidate
    row_fmiss_mref = {
        "CHROM": "chr1",
        "POS": 100,
        "proband": "0/1",
        "father": "./.",
        "mother": "0/0",
    }
    pats = deduce_patterns_for_sample("proband", row_fmiss_mref, trio_pedigree, sample_list)
    assert "de_novo_candidate" in pats
    assert "de_novo" not in pats

    # 4. 1/1 + ./. -> rejected
    row_fhom_mmiss = {
        "CHROM": "chr1",
        "POS": 100,
        "proband": "0/1",
        "father": "1/1",
        "mother": "./.",
    }
    pats = deduce_patterns_for_sample("proband", row_fhom_mmiss, trio_pedigree, sample_list)
    assert "de_novo" not in pats
    assert "de_novo_candidate" not in pats

    # 5. ./. + 1/1 -> rejected
    row_fmiss_mhom = {
        "CHROM": "chr1",
        "POS": 100,
        "proband": "0/1",
        "father": "./.",
        "mother": "1/1",
    }
    pats = deduce_patterns_for_sample("proband", row_fmiss_mhom, trio_pedigree, sample_list)
    assert "de_novo" not in pats
    assert "de_novo_candidate" not in pats

    # 6. ./. + 0/1 (mother carrier, father missing) -> rejected
    row_fmiss_mhet = {
        "CHROM": "chr1",
        "POS": 100,
        "proband": "0/1",
        "father": "./.",
        "mother": "0/1",
    }
    pats = deduce_patterns_for_sample("proband", row_fmiss_mhet, trio_pedigree, sample_list)
    assert "de_novo" not in pats
    assert "de_novo_candidate" not in pats

    # 7. ./. + ./. -> de_novo_candidate (parental data incomplete, neither carries ALT)
    row_both_miss = {
        "CHROM": "chr1",
        "POS": 100,
        "proband": "0/1",
        "father": "./.",
        "mother": "./.",
    }
    pats = deduce_patterns_for_sample("proband", row_both_miss, trio_pedigree, sample_list)
    assert "de_novo" not in pats
    assert "de_novo_candidate" in pats


@pytest.mark.unit
def test_autosomal_patterns_chrom_aliases_and_controls(trio_pedigree):
    """Verify aliases (M, CHRMT, 24, 26) are excluded and autosomes (1, 22, XY-PAR 25) included."""
    sample_list = ["proband", "father", "mother"]

    # Aliases for sex / mitochondrial chromosomes
    excluded_chroms = [
        "M",
        "chrM",
        "MT",
        "chrMT",
        "CHRMT",
        "24",
        "chr24",
        "26",
        "chr26",
        "X",
        "Y",
        " chrX ",
    ]
    for chrom in excluded_chroms:
        row = {"CHROM": chrom, "POS": 100, "proband": "0/1", "father": "0/1", "mother": "0/0"}
        pats = deduce_patterns_for_sample("proband", row, trio_pedigree, sample_list)
        for p in pats:
            assert not p.startswith("autosomal_"), f"Unexpected {p} on {chrom}"

    # Autosomal controls must receive autosomal_dominant
    autosomal_controls = ["chr1", "1", "chr2", "2", "chr22", "22", "25", "chr25"]
    for chrom in autosomal_controls:
        row = {"CHROM": chrom, "POS": 100, "proband": "0/1", "father": "0/1", "mother": "0/0"}
        df = pd.DataFrame([row])
        vec_pats = vectorized_deduce_patterns(df, trio_pedigree, sample_list)
        assert any(p.startswith("autosomal_") for p in vec_pats[0]), (
            f"Expected autosomal pattern on control {chrom}"
        )


@pytest.mark.unit
def test_mixed_compound_het_partner_types_and_patterns(trio_pedigree):
    """Verify variant with both confirmed trans and ambiguous partners
    preserves per-partner types.
    """
    import json

    from variantcentrifuge.inheritance.analyzer import analyze_inheritance

    sample_list = ["proband", "father", "mother"]

    # Var A: father-only (100)
    # Var B: mother-only (200) -> A-B is confirmed trans
    # Var C: both parents (300) -> A-C is ambiguous
    df = pd.DataFrame(
        [
            {
                "CHROM": "chr1",
                "POS": "100",
                "REF": "A",
                "ALT": "T",
                "GENE": "TEST_GENE",
                "proband": "0/1",
                "father": "0/1",
                "mother": "0/0",
            },
            {
                "CHROM": "chr1",
                "POS": "200",
                "REF": "C",
                "ALT": "G",
                "GENE": "TEST_GENE",
                "proband": "0/1",
                "father": "0/0",
                "mother": "0/1",
            },
            {
                "CHROM": "chr1",
                "POS": "300",
                "REF": "G",
                "ALT": "A",
                "GENE": "TEST_GENE",
                "proband": "0/1",
                "father": "0/1",
                "mother": "0/1",
            },
        ]
    )

    res_df = analyze_inheritance(df, trio_pedigree, sample_list)
    assert len(res_df) == 3

    # Var A details
    row_a = res_df.iloc[0]
    details_a = json.loads(row_a["Inheritance_Details"])
    assert details_a["primary_pattern"] == "compound_heterozygous"
    assert "compound_heterozygous" in details_a["all_patterns"]
    assert "compound_heterozygous_possible" in details_a["all_patterns"]

    proband_info = next(s for s in details_a["samples_with_pattern"] if s["sample_id"] == "proband")
    assert "compound_het_partner_types" in proband_info
    assert proband_info["compound_het_partner_types"]["chr1:200:C>G"] == "trans"
    assert proband_info["compound_het_partner_types"]["chr1:300:G>A"] == "ambiguous"


@pytest.mark.unit
def test_x_and_mitochondrial_case_insensitive_chrom_normalization(trio_pedigree):
    """Verify case-insensitive chromosome strings like 'ChrX' and 'chrm'
    are properly recognized in both scalar and vectorized deduction and in full
    analyze_inheritance.
    """
    import json

    from variantcentrifuge.inheritance.analyzer import analyze_inheritance

    sample_list = ["proband", "father", "mother"]

    # Proband is male affected, mother is carrier (0/1), father is ref (0/0)
    # 1. ChrX -> x_linked_recessive
    row_x = {
        "CHROM": "ChrX",
        "POS": "100",
        "REF": "A",
        "ALT": "T",
        "GENE": "X_GENE",
        "proband": "0/1",
        "father": "0/0",
        "mother": "0/1",
    }
    pats_scalar_x = deduce_patterns_for_sample("proband", row_x, trio_pedigree, sample_list)
    assert "x_linked_recessive" in pats_scalar_x

    # 2. chrm -> mitochondrial
    row_mt = {
        "CHROM": "chrm",
        "POS": "200",
        "REF": "C",
        "ALT": "G",
        "GENE": "MT_GENE",
        "proband": "0/1",
        "father": "0/0",
        "mother": "0/1",
    }
    pats_scalar_mt = deduce_patterns_for_sample("proband", row_mt, trio_pedigree, sample_list)
    assert "mitochondrial" in pats_scalar_mt

    # Full analyze_inheritance (vectorized pass)
    df = pd.DataFrame([row_x, row_mt])
    res_df = analyze_inheritance(df, trio_pedigree, sample_list)
    assert len(res_df) == 2

    # Check ChrX result
    details_x = json.loads(res_df.iloc[0]["Inheritance_Details"])
    assert details_x["primary_pattern"] == "x_linked_recessive"
    assert "x_linked_recessive" in details_x["all_patterns"]

    # Check chrm result
    details_mt = json.loads(res_df.iloc[1]["Inheritance_Details"])
    assert details_mt["primary_pattern"] == "mitochondrial"
    assert "mitochondrial" in details_mt["all_patterns"]


@pytest.mark.unit
@pytest.mark.parametrize(
    "chrom,pos,expect_ad",
    [
        ("chr7", 600_000, True),
        ("chrX", 600_000, True),  # PAR1
        ("X", 155_800_000, True),  # PAR2
        ("chrX", 2_781_480, False),  # first base after PAR1
    ],
)
def test_par_variants_get_autosomal_labels(chrom, pos, expect_ad):
    """X PAR1/PAR2 are inherited like autosomes in scalar and vectorized deduction."""
    ped = {
        "son": {"sample_id": "son", "father_id": "dad", "mother_id": "mom", "sex": "1",
                "affected_status": "2"},
        "dad": {"sample_id": "dad", "father_id": "0", "mother_id": "0", "sex": "1",
                "affected_status": "2"},
        "mom": {"sample_id": "mom", "father_id": "0", "mother_id": "0", "sex": "2",
                "affected_status": "1"},
    }  # fmt: skip
    samples = ["son", "dad", "mom"]
    row = {"CHROM": chrom, "POS": pos, "son": "0/1", "dad": "0/1", "mom": "0/0"}

    scalar = deduce_patterns_for_sample("son", row, ped, samples)
    vec = vectorized_deduce_patterns(pd.DataFrame([row]), ped, samples)[0]
    assert ("autosomal_dominant" in scalar) is expect_ad
    assert ("autosomal_dominant" in vec) is expect_ad
    if expect_ad:
        assert not any(p.startswith("x_linked") for p in scalar + vec)


@pytest.mark.unit
def test_chrom_25_is_not_mitochondrial(trio_pedigree):
    """PLINK 25 is the XY PAR, 26 is MT."""
    samples = ["proband", "father", "mother"]
    for chrom, expect_mt in [("25", False), ("chr25", False), ("26", True), ("chrM", True)]:
        row = {"CHROM": chrom, "POS": 100, "proband": "0/1", "father": "0/0", "mother": "0/1"}
        scalar = deduce_patterns_for_sample("proband", row, trio_pedigree, samples)
        vec = vectorized_deduce_patterns(pd.DataFrame([row]), trio_pedigree, samples)[0]
        assert ("mitochondrial" in scalar) is expect_mt, chrom
        assert ("mitochondrial" in vec) is expect_mt, chrom


@pytest.mark.unit
@pytest.mark.parametrize("mother_gt", ["0/0", "0/1", "1/1"])
@pytest.mark.parametrize("child_sex", ["1", "2"])
def test_x_linked_scalar_vectorized_parity_single_known_mother(mother_gt, child_sex):
    """Non-PAR X, father unknown ("0"): scalar and vectorized agree.

    Regression: vectorized XLD skipped its own "x_linked_dominant_possible" because another
    sample had already put "x_linked_dominant" in the shared per-variant list.
    """
    from variantcentrifuge.inheritance.deducer import deduce_patterns_for_variant

    ped = {
        "kid": {"sample_id": "kid", "father_id": "0", "mother_id": "mom", "sex": child_sex,
                "affected_status": "2"},
        "mom": {"sample_id": "mom", "father_id": "0", "mother_id": "0", "sex": "2",
                "affected_status": "2"},
    }  # fmt: skip
    samples = ["kid", "mom"]
    row = {"CHROM": "chrX", "POS": 50_000_000, "kid": "0/1", "mom": mother_gt}
    scalar = deduce_patterns_for_variant(row, ped, samples)
    vec = vectorized_deduce_patterns(pd.DataFrame([row]), ped, samples)[0]
    assert sorted(scalar) == sorted(vec)
