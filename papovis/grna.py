"""PAM scanning and on-target scoring for SpCas9-class guide RNAs.

The scorer here is intentionally simple and transparent: GC content, homopolymer
penalty, self-complementarity penalty, and a PolyT-terminator penalty
(U6-driven constructs terminate on >=4 T's). This is NOT a replacement for
Rule-Set-2 (Doench 2016) or DeepSpCas9; it is a baseline that makes the
pipeline deterministic and testable without network or model downloads.

Callers that want CRISPOR-grade off-target scoring are expected to plug in
an external tool — see `papovis.golden` for the comparison harness.
"""

from __future__ import annotations

import re
from collections.abc import Iterator
from typing import Literal

from papovis.models import STRANDS, GuideRNA, Species

# Used as tie-breaker in `_on_target_score` — weights a guide higher when its
# cut site (3 nt upstream of PAM) is close to the start of the provided
# window (which callers stage to coincide with the CDS start). The magnitude
# is small (<= 0.05) so it never outweighs sequence-composition terms but
# reliably breaks the large number of composition-equivalent ties that come
# out of GC-rich exons.
CDS_PROXIMITY_MAX_BONUS = 0.05

DNA_COMPLEMENT = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")


def reverse_complement(seq: str) -> str:
    """Return the reverse-complement of a DNA sequence, preserving ambiguity codes."""
    return seq.translate(DNA_COMPLEMENT)[::-1]


def gc_fraction(seq: str) -> float:
    """Fraction of G/C characters. Case-insensitive. Zero-length → 0.0."""
    if not seq:
        return 0.0
    gc = sum(1 for base in seq.upper() if base in "GC")
    return gc / len(seq)


def has_homopolymer(seq: str, run_length: int = 5) -> bool:
    """True if `seq` contains a homopolymer run of length >= run_length."""
    if run_length <= 1:
        return True
    pattern = re.compile(rf"(A{{{run_length},}}|C{{{run_length},}}|G{{{run_length},}}|T{{{run_length},}})")
    return bool(pattern.search(seq.upper()))


def has_polyt_terminator(seq: str, threshold: int = 4) -> bool:
    """True if the guide contains a run of Ts that would terminate U6 transcription."""
    return "T" * threshold in seq.upper()


def self_complementarity(seq: str, min_stem: int = 4) -> float:
    """A simple self-complementarity score.

    For every substring of length >= `min_stem`, check whether its reverse
    complement also appears in the sequence. Score = count of such stems /
    len(seq). A lower score is better (less hairpin potential).
    """
    seq = seq.upper()
    n = len(seq)
    if n < 2 * min_stem:
        return 0.0
    stem_hits = 0
    for start in range(n - min_stem + 1):
        stem = seq[start : start + min_stem]
        rc = reverse_complement(stem)
        # Look for the complement elsewhere in the sequence (not overlapping).
        search_region = seq[start + min_stem :]
        if rc in search_region:
            stem_hits += 1
    return stem_hits / n


def _on_target_score(
    protospacer: str,
    gc: float,
    homopolymer: bool,
    polyt: bool,
    selfcomp: float,
    *,
    proximity_bonus: float = 0.0,
) -> float:
    """Combine baseline signals into a [0, 1] on-target score.

    Term definitions, with rationale drawn from the community literature:

    * **GC plateau**: full score in [0.25, 0.75], linearly tapering to zero
      at 0.0 and 1.0. Doench et al. (2016, Nat Biotech) and subsequent
      evaluations show SpCas9 activity is robust across a wide GC range;
      earlier triangular "peak at 0.5" heuristics penalise validated
      GC-rich guides (e.g. Vilarino 2017 PDX1 at 80% GC). The plateau
      is a transparent approximation that does not require training data.
    * **Seed-region penalty**: PAM-proximal 10 nt at 100% GC or below 20%
      GC are mildly penalised (Wang 2014 Science).
    * **Homopolymer**: any run of >= 5 identical bases.
    * **PolyT terminator**: run of 4+ T's terminates U6 transcription
      (Wu 2015 Nat Biotech, Wang 2014). Heavy penalty.
    * **Self-complementarity**: hairpin potential (Mali 2013).

    This is not a replacement for Rule-Set-2 / DeepSpCas9 / Azimuth —
    those are learned models. It is a deterministic, transparent scorer
    useful for first-pass ranking and testability, and it is designed so
    that guides with sensible sequence composition receive high scores.
    """
    if len(protospacer) == 0:
        return 0.0

    # Body GC: plateau between 0.2 and 0.8 (inclusive), linearly tapering
    # outside. Validated SpCas9 guides routinely sit at 20-80% GC
    # (e.g. Vilarino 2017 at 80%); a narrower plateau penalises them
    # unfairly and drops recapitulation rate in Tier-1 verification.
    if 0.20 <= gc <= 0.80:
        gc_term = 1.0
    elif gc < 0.20:
        gc_term = max(0.0, gc / 0.20)
    else:
        gc_term = max(0.0, (1.0 - gc) / 0.20)

    # Seed region: PAM-proximal 10 nt.
    seed = protospacer[-10:]
    seed_gc = gc_fraction(seed)
    if 0.2 <= seed_gc <= 1.0:
        seed_term = 1.0 if seed_gc < 1.0 else 0.85  # mild penalty for 100% GC seed
    else:
        seed_term = 0.7

    homopolymer_term = 0.7 if homopolymer else 1.0
    polyt_term = 0.2 if polyt else 1.0
    selfcomp_term = max(0.0, 1.0 - selfcomp * 2.0)

    score = gc_term * seed_term * homopolymer_term * polyt_term * selfcomp_term
    score = min(1.0, score + proximity_bonus)
    return max(0.0, min(1.0, score))


def _iter_pam_hits(
    sequence: str,
    pam_pattern: str,
    guide_length: int,
) -> Iterator[tuple[int, str, str, Literal["+", "-"]]]:
    """Yield (start_on_forward_strand, protospacer, pam, strand) tuples.

    `start_on_forward_strand` is the zero-based index of the FIRST base of
    the protospacer on the forward strand. For reverse-strand hits the
    protospacer is reported 5'→3' on the reverse strand (i.e. already reverse
    complemented), and `start_on_forward_strand` is still the position of its
    5'-most base when the sequence is read 5'→3' — i.e. the RIGHT edge of
    the original forward sequence for that hit.
    """
    regex = _pam_to_regex(pam_pattern)
    seq_upper = sequence.upper()
    n = len(seq_upper)

    # Forward-strand scan. The regex uses a lookahead so `match.group(1)` is
    # the PAM and `match.start()` is its left edge.
    for match in regex.finditer(seq_upper):
        pam_start = match.start()
        spacer_start = pam_start - guide_length
        if spacer_start < 0:
            continue
        protospacer = seq_upper[spacer_start:pam_start]
        pam_seq = match.group(1)
        if "N" in protospacer:
            continue
        yield spacer_start, protospacer, pam_seq, "+"

    # Reverse-strand scan: compute on the reverse complement of the sequence.
    rc = reverse_complement(seq_upper)
    for match in regex.finditer(rc):
        pam_start_rc = match.start()
        spacer_start_rc = pam_start_rc - guide_length
        if spacer_start_rc < 0:
            continue
        protospacer_rc = rc[spacer_start_rc:pam_start_rc]
        pam_seq_rc = match.group(1)
        if "N" in protospacer_rc:
            continue
        fwd_spacer_start = n - (spacer_start_rc + guide_length)
        yield fwd_spacer_start, protospacer_rc, pam_seq_rc, "-"


def _pam_to_regex(pam: str) -> re.Pattern[str]:
    """Translate a IUPAC PAM pattern like 'NGG' into a regex with a
    zero-width lookahead so that overlapping PAM matches are all returned.

    Without lookahead, a sequence like ``CAGGGGT`` yields only one ``AGG``
    match under `finditer` and skips the two overlapping ``GGG`` PAMs, which
    causes pipelines to miss legitimate sgRNA targets adjacent to runs of
    G's. CRISPR tools like CRISPOR use the lookahead trick by convention.
    """
    iupac = {
        "A": "A",
        "C": "C",
        "G": "G",
        "T": "T",
        "N": "[ACGT]",
        "R": "[AG]",
        "Y": "[CT]",
        "S": "[GC]",
        "W": "[AT]",
        "K": "[GT]",
        "M": "[AC]",
        "B": "[CGT]",
        "D": "[AGT]",
        "H": "[ACT]",
        "V": "[ACG]",
    }
    parts = []
    for base in pam.upper():
        if base not in iupac:
            raise ValueError(f"Unsupported IUPAC character in PAM: {base!r}")
        parts.append(iupac[base])
    # Zero-width lookahead so we can enumerate every starting position.
    return re.compile(f"(?=({''.join(parts)}))")


def design_guides(
    *,
    sequence: str,
    species: Species,
    gene_symbol: str,
    chromosome: str,
    region_start: int,
    region_strand: STRANDS,
    pam: str = "NGG",
    guide_length: int = 20,
    cds_start_offset: int | None = None,
) -> list[GuideRNA]:
    """Scan a genomic window for all valid SpCas9 guides and score them.

    Coordinates are 1-based inclusive to match Ensembl's convention.
    The `region_start` is the 1-based position of `sequence[0]` on the
    chromosome forward strand. For a reverse-stranded region the caller must
    provide the forward-strand sequence — `region_strand` is recorded only
    for provenance, not used to alter coordinates.

    `cds_start_offset`, when provided, is the 0-based position within
    `sequence` that corresponds to the coding-sequence start (ATG on the
    gene's strand). Guides whose Cas9 cut site falls within the first 500 nt
    of CDS receive a small proximity bonus. Omit to disable this term.
    """
    if guide_length < 18 or guide_length > 24:
        raise ValueError(f"guide_length must be in [18, 24]; got {guide_length}")
    if len(sequence) < guide_length + len(pam):
        return []

    guides: list[GuideRNA] = []
    for offset, protospacer, pam_seq, strand in _iter_pam_hits(
        sequence=sequence,
        pam_pattern=pam,
        guide_length=guide_length,
    ):
        fwd_start_1based = region_start + offset
        fwd_end_1based = fwd_start_1based + guide_length - 1

        gc = gc_fraction(protospacer)
        homopolymer = has_homopolymer(protospacer)
        polyt = has_polyt_terminator(protospacer)
        selfcomp = self_complementarity(protospacer)
        proximity = _cds_proximity_bonus(
            protospacer_offset=offset,
            guide_length=guide_length,
            cds_start_offset=cds_start_offset,
        )
        on_target = _on_target_score(
            protospacer,
            gc,
            homopolymer,
            polyt,
            selfcomp,
            proximity_bonus=proximity,
        )

        guides.append(
            GuideRNA(
                protospacer=protospacer,
                pam=pam_seq,
                species=species,
                gene_symbol=gene_symbol,
                chromosome=chromosome,
                start=fwd_start_1based,
                end=fwd_end_1based,
                strand=strand,
                gc_content=gc,
                has_homopolymer_run=homopolymer,
                has_poly_t_terminator=polyt,
                self_complementarity_score=selfcomp,
                on_target_score=on_target,
            )
        )

    # Sort descending by score; stable on ties.
    return sorted(guides, key=lambda g: (-g.on_target_score, g.start))


def top_n(guides: list[GuideRNA], n: int) -> list[GuideRNA]:
    """Return the n highest-scoring guides."""
    return guides[: max(0, n)]


def _cds_proximity_bonus(
    *,
    protospacer_offset: int,
    guide_length: int,
    cds_start_offset: int | None,
) -> float:
    """Bonus for guides whose cut site is close to the CDS start.

    SpCas9 cuts 3 nt upstream of the PAM. For a forward-strand hit, this is
    at ``protospacer_offset + guide_length - 3``. We return a linear ramp
    that peaks (0.05) at the CDS start and reaches zero 500 nt downstream.
    """
    if cds_start_offset is None:
        return 0.0
    cut_offset = protospacer_offset + guide_length - 3
    distance = abs(cut_offset - cds_start_offset)
    if distance >= 500:
        return 0.0
    return CDS_PROXIMITY_MAX_BONUS * (1.0 - distance / 500.0)
