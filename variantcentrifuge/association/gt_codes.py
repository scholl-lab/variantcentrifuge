"""Run-level genotype-code matrix: parse genotypes once, slice them per gene (#219).

The association stage parses the packed ``GT`` column (``Sample(0/1);...`` or positional
``0/1,./.,...``), or the per-sample ``GEN[i].GT`` columns when no packed ``GT`` exists, ONCE
into an integer matrix ``codes`` (variants x analysed samples). Each code indexes ``vocab``,
the distinct genotype strings, so per-string semantics are classified once per vocabulary
entry and every per-cell property is a table lookup:

- the Fisher/collapsing counts (:func:`aggregate_gene_counts`), and
- the regression dosage/ploidy/kind arrays (:meth:`GenotypeCodes.gene_arrays`) consumed by
  :func:`~variantcentrifuge.association.genotype_matrix.build_genotype_matrix`.

Packed rows list only non-reference samples; unlisted samples are ``0/0``.
"""

from __future__ import annotations

import logging
from collections.abc import Iterator, Sequence
from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd
import pyarrow as pa  # type: ignore[import-untyped]
import pyarrow.compute as pc  # type: ignore[import-untyped]

from ..genotype_utils import is_missing
from ..gt_reconstruction import gt_columns_match_samples
from .aggregation import _gt_to_dosage
from .engine_workers import fork_available, fork_map
from .genotype_matrix import _classify_gt_vocabulary, _factorize_gt_block, _multiallelic_error

logger = logging.getLogger("variantcentrifuge")

# "Sample(GT[:extra...])": name before "(", GT = first token inside, like
# phenotype_helpers.extract_sample_and_genotype (strip, split on [:,;\s]).
_ENTRY = r"^\s*(?P<name>[^(]*?)\s*\(\s*(?P<gt>[^:,;\s)]*)[^)]*\)"
_EMPTY_ROWS = {"", "NA", "nan"}
_BATCH_CHARS = 64_000_000  # GT text per Arrow batch: transient arrays are ~7x the text


class GenotypeCodes:
    """Genotype codes for all variants x analysed samples plus per-code lookup tables.

    ``codes`` is int8 (int16 only if more than 127 distinct genotype strings occur).
    Builders keep a reference to one shared instance and slice rows at call time.
    """

    def __init__(self, codes: np.ndarray, vocab: list[str], labels: list[str]) -> None:
        self.codes = codes
        self.vocab = vocab
        self.labels = labels  # per sample: source column (or sample ID) for error messages
        s = pd.Series(vocab, dtype=object)
        self.dosage, self.ploidy, self.kind, self.multi = _classify_gt_vocabulary(s)
        stripped = [v.strip() for v in vocab]
        self.missing = np.array([is_missing(v) for v in stripped], dtype=bool)
        self.dose = np.array([_gt_to_dosage(v) for v in vocab], dtype=np.int8)
        # Ploidy vote for the collapsing counts: 1 haploid, 2 diploid, 0 none (empty cell).
        self.vote = np.array(
            [
                (0 if v == "" else 1 if v == "." else 2)
                if m
                else (2 if "/" in v or "|" in v else 1)
                for v, m in zip(stripped, self.missing, strict=True)
            ],
            dtype=np.int8,
        )

    def gene_arrays(
        self, rows: np.ndarray, gene_df: pd.DataFrame
    ) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """Return (dosage, ploidy, kind) for *rows*, exactly as the column parser does."""
        sub = self.codes[rows]
        multi = self.multi[sub] if self.multi.any() else None
        if multi is not None and multi.any():  # first offending sample, then row (column-major)
            j, r = divmod(int(np.flatnonzero(multi.T)[0]), len(rows))
            raise _multiallelic_error(gene_df, r, self.labels[j], self.vocab[sub[r, j]])
        return self.dosage[sub], self.ploidy[sub], self.kind[sub]


class _Vocab:
    """Growing string -> code table; codes start as int8 and widen if needed."""

    def __init__(self) -> None:
        self.index: dict[str, int] = {"0/0": 0}
        self.items: list[str] = ["0/0"]

    def lut(self, strings: Sequence[str]) -> np.ndarray:
        for s in strings:
            if s not in self.index:
                self.index[s] = len(self.items)
                self.items.append(s)
        return np.array([self.index[s] for s in strings], dtype=np.int16)


def _store(codes: np.ndarray, vocab: _Vocab) -> np.ndarray:
    """Widen int8 codes to int16 once the vocabulary outgrows int8."""
    if codes.dtype == np.int8 and len(vocab.items) > 127:
        return codes.astype(np.int16)
    return codes


def _prefault(codes: np.ndarray, workers: int) -> None:
    """Touch the pages of a fresh ``np.zeros`` matrix from ``workers`` threads.

    A run-level code matrix is gigabytes; first-touch page faults from the merge's
    scattered writes cost seconds in one thread, and threads fault pages in parallel.
    Call only once the parse workers are forked: pages populated before a fork would be
    copied on write while the children live.
    """
    step = max(1, -(-codes.shape[0] // (workers * 4)))
    with ThreadPoolExecutor(workers) as pool:  # ndarray.fill releases the GIL
        list(pool.map(lambda lo: codes[lo : lo + step].fill(0), range(0, codes.shape[0], step)))


def parse_genotype_codes(
    df: pd.DataFrame,
    samples: list[str],
    vcf_samples: list[str] | None,
    gt_columns: list[str],
    workers: int = 1,
) -> GenotypeCodes:
    """Parse per-sample GT columns (preferred) or the packed ``GT`` for *samples* once.

    Columns are used only when they are exactly ``GEN[0..n-1].GT`` for the n VCF samples,
    in order (``gt_columns[i]`` holds ``vcf_samples[i]``). They win over packed GT because
    every sample ID survives in them, while IDs containing ``(``/``;`` or edge whitespace
    cannot round-trip through ``Name(GT);...``. ``workers`` > 1 parses packed GT in
    forked processes; the result is the same for any worker count.
    """
    vocab = _Vocab()
    codes = np.zeros((len(df), len(samples)), dtype=np.int8)
    columns_ok = gt_columns_match_samples(gt_columns, vcf_samples)
    if not columns_ok and "GT" not in df.columns:
        raise ValueError(
            f"{len(gt_columns)} per-sample GT columns do not map 1:1 (GEN[0..n-1].GT, in "
            f"order) to the {len(vcf_samples or [])} VCF samples and there is no packed GT "
            "column; cannot assign genotypes to samples"
        )
    if not columns_ok and gt_columns:
        logger.warning(
            "%d per-sample GT columns do not map 1:1 to the %d VCF samples; using packed GT",
            len(gt_columns),
            len(vcf_samples or []),
        )
    if not columns_ok:
        ids = vcf_samples or samples
        bad = sum(1 for s in ids if "(" in s or ";" in s or s != s.strip())
        if bad:
            logger.warning(
                "%d of %d VCF sample IDs cannot be represented in packed GT (they contain "
                "'(' or ';' or leading/trailing spaces); their genotypes read as reference",
                bad,
                len(ids),
            )
        codes = _parse_packed(df["GT"], samples, vcf_samples, codes, vocab, workers)
        labels = list(samples)
    else:
        col_of = dict(zip(vcf_samples or [], gt_columns, strict=True))
        labels = [col_of[s] for s in samples]
        codes = _parse_columns(df, labels, codes, vocab)
    return GenotypeCodes(codes, vocab.items, labels)


def _parse_columns(
    df: pd.DataFrame, columns: list[str], codes: np.ndarray, vocab: _Vocab
) -> np.ndarray:
    """Encode per-sample GT columns; cells normalised like ``fillna("").astype(str)``."""
    for j, col in enumerate(columns):
        if isinstance(df[col].dtype, pd.CategoricalDtype):
            series = df[col]
            lut = vocab.lut([str(c) for c in series.cat.categories] + [""])
            cell = series.cat.codes.to_numpy()  # -1 (NA) picks the trailing ""
        else:
            cell, strings = _factorize_gt_block(df[[col]])
            lut = vocab.lut(strings)
        codes = _store(codes, vocab)
        codes[:, j] = lut[cell]
    return codes


def _char_batches(text: np.ndarray) -> Iterator[tuple[int, int]]:
    """Row ranges of about ``_BATCH_CHARS`` characters of GT text (at least one row)."""
    ends = np.cumsum(np.fromiter(map(len, text), dtype=np.int64, count=len(text)))
    lo = 0
    while lo < len(text):
        start = int(ends[lo - 1]) if lo else 0
        hi = max(lo + 1, int(np.searchsorted(ends, start + _BATCH_CHARS, side="right")))
        yield lo, hi
        lo = hi


def _parse_packed_rows(
    state: tuple[np.ndarray, pa.Array], lo: int, hi: int
) -> tuple[np.ndarray, np.ndarray, np.ndarray, list[str], np.ndarray]:
    """Parse packed GT rows text[lo:hi] (a ``fork_map`` task; also the sequential path).

    Arrow work runs in batches of about ``_BATCH_CHARS`` characters. Returns (rows naming
    a sample, then for entries of wanted samples: absolute row, sample column, the
    distinct GT strings in order of appearance and each entry's index into them).
    """
    text, wanted = state
    local = _Vocab()
    parts: list[tuple[np.ndarray, ...]] = []
    for b_lo, b_hi in _char_batches(text[lo:hi]):
        chunk = pa.array(text[lo + b_lo : lo + b_hi], type=pa.large_string())
        named = pc.match_substring(chunk, "(").to_numpy(False)
        lists = pc.split_pattern(chunk, ";")
        rows = pc.list_parent_indices(lists).to_numpy() + lo + b_lo
        entries = pc.extract_regex(pc.list_flatten(lists), _ENTRY)
        cols = pc.fill_null(pc.index_in(entries.field("name"), value_set=wanted), -1)
        cols = cols.to_numpy(zero_copy_only=False)
        keep = cols >= 0
        enc = pc.dictionary_encode(pc.filter(entries.field("gt"), pa.array(keep)))
        enc = enc.combine_chunks() if isinstance(enc, pa.ChunkedArray) else enc
        lut = local.lut(enc.dictionary.to_pylist())
        parts.append((named, rows[keep], cols[keep], lut[enc.indices.to_numpy()]))
    named, rows, cols, idx = (np.concatenate(p) for p in zip(*parts, strict=True))
    return named, rows, cols, local.items, idx


def _parse_packed(
    gt: pd.Series,
    samples: list[str],
    vcf_samples: list[str] | None,
    codes: np.ndarray,
    vocab: _Vocab,
    workers: int = 1,
) -> np.ndarray:
    """Encode a packed GT column for *samples* (unlisted samples stay ``0/0``).

    With ``workers`` > 1 the rows are parsed in forked processes (``fork_map``);
    only the wanted entries come back, merged here in row order.
    """
    if isinstance(gt.dtype, pd.CategoricalDtype):
        gt = gt.astype(object)
    text = gt.fillna("").astype(str).to_numpy(dtype=object)
    state = (text, pa.array(samples, type=pa.large_string()))
    named = np.zeros(len(text), dtype=bool)
    forked = workers > 1 and len(text) > 1 and fork_available()
    if forked:
        parts = fork_map(_parse_packed_rows, state, len(text), workers)
    else:
        # One batch per call, merged as it comes: no run-wide list of entries.
        parts = ((lo, hi, _parse_packed_rows(state, lo, hi)) for lo, hi in _char_batches(text))
    for lo, hi, (row_named, rows, cols, gts, idx) in parts:
        if forked and lo == 0:  # first result: every worker is forked by now
            _prefault(codes, workers)
        named[lo:hi] = row_named
        lut = vocab.lut(gts)
        codes = _store(codes, vocab)
        codes[rows, cols] = lut[idx]
    positional = ~named & ~np.isin(text, list(_EMPTY_ROWS))
    if positional.any() and vcf_samples:
        # Positional rows list every VCF sample in order; absent trailing calls are ./.
        pos_of = {s: i for i, s in enumerate(vcf_samples)}
        split = pd.Series(text[positional]).str.split(",", expand=True)
        for j, s in enumerate(samples):
            i = pos_of.get(s)
            values = split[i] if i is not None and i in split.columns else None
            strings = ["./."] * int(positional.sum()) if values is None else values
            cell, uniq = pd.factorize(pd.Series(strings, dtype=object).fillna("./."))
            lut = vocab.lut(list(uniq))
            codes = _store(codes, vocab)
            codes[positional, j] = lut[cell]
    return codes


def aggregate_gene_counts(
    genotypes: GenotypeCodes,
    gene_rows: dict[str, np.ndarray],
    case_mask: np.ndarray,
    control_mask: np.ndarray,
    missing_as_ref: bool = False,
    workers: int = 1,
) -> list[dict]:
    """Per-gene case/control collapsing counts (Fisher / burden tables) from the codes.

    Complete-case by default: a sample with any missing call in the gene is not counted;
    ``missing_as_ref`` counts missing calls as reference (#200). A called sample is
    haploid when its calls vote haploid and none votes diploid. ``workers`` > 1 counts
    gene ranges in forked processes; the result is the same for any worker count.
    """
    state = (genotypes, list(gene_rows.items()), case_mask, control_mask, missing_as_ref)
    n = len(state[1])
    if workers > 1 and n > 1 and fork_available():
        return [row for _lo, _hi, part in fork_map(_count_genes, state, n, workers) for row in part]
    return _count_genes(state, 0, n)


def _count_genes(state: tuple, lo: int, hi: int) -> list[dict]:
    """Counts for genes[lo:hi] (a ``fork_map`` task; also the sequential loop)."""
    genotypes, items, case_mask, control_mask, missing_as_ref = state
    out = []
    n_cases, n_controls = int(case_mask.sum()), int(control_mask.sum())
    for gene, rows in items[lo:hi]:
        sub = genotypes.codes[rows]
        called = (
            np.ones(sub.shape[1], dtype=bool)
            if missing_as_ref
            else ~genotypes.missing[sub].any(axis=0)
        )
        dose = genotypes.dose[sub].max(axis=0).astype(np.int64)
        chrom = np.where(genotypes.vote[sub].max(axis=0) == 1, 1, 2)
        carrier = called & (dose > 0)
        counts = {}
        for prefix, mask in (("proband", case_mask), ("control", control_mask)):
            m_called = called & mask
            counts[f"{prefix}_called_count"] = int(m_called.sum())
            counts[f"{prefix}_called_chromosomes"] = int(chrom[m_called].sum())
            counts[f"{prefix}_carrier_count"] = int((carrier & mask).sum())
            counts[f"{prefix}_allele_count"] = int(dose[carrier & mask].sum())
        out.append(
            {
                "GENE": gene,
                "proband_count": n_cases,
                "control_count": n_controls,
                "proband_called_count": counts["proband_called_count"],
                "control_called_count": counts["control_called_count"],
                "proband_called_chromosomes": counts["proband_called_chromosomes"],
                "control_called_chromosomes": counts["control_called_chromosomes"],
                "proband_carrier_count": counts["proband_carrier_count"],
                "control_carrier_count": counts["control_carrier_count"],
                "proband_allele_count": counts["proband_allele_count"],
                "control_allele_count": counts["control_allele_count"],
                "n_qualifying_variants": len(rows),
            }
        )
    return out
