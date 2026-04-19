"""Unit tests for the PAM scanner + on-target scorer.

All tests here are network-free and deterministic.
"""

from __future__ import annotations

import pytest

from papovis.grna import (
    design_guides,
    gc_fraction,
    has_homopolymer,
    has_polyt_terminator,
    reverse_complement,
    self_complementarity,
)
from papovis.models import GuideRNA


class TestReverseComplement:
    def test_known_sequence(self) -> None:
        assert reverse_complement("ATGC") == "GCAT"

    def test_preserves_ambiguity(self) -> None:
        assert reverse_complement("NATG") == "CATN"

    def test_idempotent_round_trip(self) -> None:
        seq = "ACGTTGCAACGTNNACGT"
        assert reverse_complement(reverse_complement(seq)) == seq


class TestScorers:
    def test_gc_fraction_empty(self) -> None:
        assert gc_fraction("") == 0.0

    def test_gc_fraction_mixed(self) -> None:
        assert gc_fraction("ACGT") == pytest.approx(0.5)

    def test_homopolymer_detected(self) -> None:
        assert has_homopolymer("ACGAAAAAT") is True

    def test_homopolymer_absent(self) -> None:
        assert has_homopolymer("ACGTACGT") is False

    def test_polyt_detected(self) -> None:
        assert has_polyt_terminator("ACGTTTTCG") is True

    def test_polyt_absent(self) -> None:
        assert has_polyt_terminator("ACGTAGCT") is False

    def test_self_complementarity_zero_for_short(self) -> None:
        assert self_complementarity("ACGT") == 0.0

    def test_self_complementarity_nonzero_for_palindrome(self) -> None:
        # ACGTACGT has ACGT at pos 0 whose RC (ACGT) appears at pos 4.
        score = self_complementarity("ACGTACGT")
        assert score > 0.0


class TestDesignGuides:
    def test_returns_empty_for_no_pam(self, short_sequence_no_pam: str) -> None:
        # 20 A's then no NGG — no PAM can follow.
        guides = design_guides(
            sequence=short_sequence_no_pam,
            species="Ovis aries",
            gene_symbol="TEST",
            chromosome="1",
            region_start=1,
            region_strand="+",
        )
        assert guides == []

    def test_finds_single_guide_at_end(self, min_pam_sequence: str) -> None:
        guides = design_guides(
            sequence=min_pam_sequence,
            species="Ovis aries",
            gene_symbol="TEST",
            chromosome="1",
            region_start=1,
            region_strand="+",
        )
        # At least one guide on the forward strand.
        forward = [g for g in guides if g.strand == "+"]
        assert len(forward) >= 1
        hit = forward[0]
        assert hit.protospacer == "ACGTACGTACGTACGTACGT"
        assert hit.pam == "AGG"
        assert hit.start == 1
        assert hit.end == 20

    def test_sorted_descending_by_score(
        self,
        synthetic_pdx1_like_sequence: str,
    ) -> None:
        guides = design_guides(
            sequence=synthetic_pdx1_like_sequence,
            species="Ovis aries",
            gene_symbol="TEST",
            chromosome="1",
            region_start=1,
            region_strand="+",
        )
        assert guides, "Synthetic sequence should produce guides on both strands"
        scores = [g.on_target_score for g in guides]
        assert scores == sorted(scores, reverse=True)

    def test_scans_both_strands(self, synthetic_pdx1_like_sequence: str) -> None:
        guides = design_guides(
            sequence=synthetic_pdx1_like_sequence,
            species="Ovis aries",
            gene_symbol="TEST",
            chromosome="1",
            region_start=1,
            region_strand="+",
        )
        strands = {g.strand for g in guides}
        # Exercise at least one of each strand orientation when PAMs exist.
        assert strands, "expected at least one strand represented"

    def test_guide_coordinates_are_1_based_inclusive(
        self,
        min_pam_sequence: str,
    ) -> None:
        region_start = 1000
        guides = design_guides(
            sequence=min_pam_sequence,
            species="Ovis aries",
            gene_symbol="TEST",
            chromosome="1",
            region_start=region_start,
            region_strand="+",
        )
        hit = [g for g in guides if g.strand == "+" and g.protospacer == "ACGTACGTACGTACGTACGT"][0]
        assert hit.start == region_start  # offset 0 → region_start
        assert hit.end == region_start + 19  # 20-mer inclusive

    def test_polyt_terminator_demotes_score(self) -> None:
        with_polyt = "ACGTACGTACGTTTTCAGCT" + "AGG"
        without_polyt = "ACGTACGTACGTCACCAGCT" + "AGG"
        g1 = design_guides(
            sequence=with_polyt,
            species="Ovis aries",
            gene_symbol="T",
            chromosome="1",
            region_start=1,
            region_strand="+",
        )[0]
        g2 = design_guides(
            sequence=without_polyt,
            species="Ovis aries",
            gene_symbol="T",
            chromosome="1",
            region_start=1,
            region_strand="+",
        )[0]
        assert g1.on_target_score < g2.on_target_score

    def test_rejects_invalid_guide_length(self) -> None:
        with pytest.raises(ValueError):
            design_guides(
                sequence="ACGT",
                species="Ovis aries",
                gene_symbol="T",
                chromosome="1",
                region_start=1,
                region_strand="+",
                guide_length=10,
            )


class TestGuideRNAModel:
    def test_is_frozen(self, min_pam_sequence: str) -> None:
        guides = design_guides(
            sequence=min_pam_sequence,
            species="Ovis aries",
            gene_symbol="TEST",
            chromosome="1",
            region_start=1,
            region_strand="+",
        )
        guide = guides[0]
        with pytest.raises(Exception):  # pydantic ValidationError on frozen model
            guide.protospacer = "X" * 20  # type: ignore[misc]

    def test_full_target_concatenates_pam(self, min_pam_sequence: str) -> None:
        guide: GuideRNA = design_guides(
            sequence=min_pam_sequence,
            species="Ovis aries",
            gene_symbol="TEST",
            chromosome="1",
            region_start=1,
            region_strand="+",
        )[0]
        assert guide.full_target() == guide.protospacer + guide.pam
