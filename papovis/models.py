"""Shared immutable pydantic models used across pipeline stages."""

from __future__ import annotations

from enum import Enum
from typing import Literal

from pydantic import BaseModel, ConfigDict, Field

Species = Literal["Ovis aries", "Papio anubis", "Homo sapiens", "Sus scrofa"]

# Ensembl genome assembly identifiers pinned at release 115 (Sep 2025).
ASSEMBLY_BY_SPECIES: dict[Species, str] = {
    "Ovis aries": "ARS-UI_Ramb_v2.0",
    "Papio anubis": "Panubis1.0",
    "Homo sapiens": "GRCh38",
    "Sus scrofa": "Sscrofa11.1",
}

STRANDS = Literal["+", "-"]


class _Frozen(BaseModel):
    """Base class for immutable pydantic models."""

    model_config = ConfigDict(frozen=True, extra="forbid")


class EditDirection(str, Enum):
    """How a gene should be manipulated."""

    KO = "KO"
    OVEREXPRESS = "OVEREXPRESS"
    KNOCK_IN = "KNOCK_IN"
    ASSESS = "ASSESS"


class EditNecessity(str, Enum):
    """Whether an xeno-antigen edit is needed for the sheep→baboon axis."""

    REQUIRED = "REQUIRED"
    RECOMMENDED = "RECOMMENDED"
    NOT_REQUIRED = "NOT_REQUIRED"
    COUNTERPRODUCTIVE = "COUNTERPRODUCTIVE"
    ASSESS = "ASSESS"


class Gene(_Frozen):
    """A gene in a given species."""

    symbol: str
    species: Species
    ensembl_id: str | None = None
    human_ortholog_ensembl_id: str | None = None
    notes: str | None = None


class SequenceRegion(_Frozen):
    """A genomic interval plus its DNA sequence on the forward strand."""

    species: Species
    chromosome: str
    start: int = Field(ge=1, description="1-based inclusive start")
    end: int = Field(ge=1, description="1-based inclusive end")
    strand: STRANDS
    sequence: str

    def length(self) -> int:
        return self.end - self.start + 1


class GuideRNA(_Frozen):
    """A candidate 20 nt spacer, recorded as 5'→3' on its target strand."""

    protospacer: str = Field(min_length=18, max_length=24)
    pam: str
    species: Species
    gene_symbol: str
    # Genomic coordinates of the protospacer (not including PAM)
    chromosome: str
    start: int
    end: int
    strand: STRANDS
    gc_content: float = Field(ge=0.0, le=1.0)
    has_homopolymer_run: bool
    has_poly_t_terminator: bool
    self_complementarity_score: float = Field(ge=0.0)
    on_target_score: float = Field(ge=0.0, le=1.0)

    def full_target(self) -> str:
        """Protospacer + PAM, 5'→3'."""
        return f"{self.protospacer}{self.pam}"


class NicheEditPlan(_Frozen):
    """Stage 1 niche-vacancy editing plan for the host (sheep) embryo."""

    organ: str
    host_species: Species
    target_genes: tuple[str, ...]
    guides: tuple[GuideRNA, ...]
    references: tuple[str, ...]


class CompetenceEditPlan(_Frozen):
    """Stage 1 donor (baboon) iPSC editing plan for chimera competence."""

    donor_species: Species
    host_species: Species
    barrier_class: str
    direction: EditDirection
    gene_symbol: str
    guides: tuple[GuideRNA, ...]
    rationale: str
    references: tuple[str, ...]


class XenoEditRow(_Frozen):
    """One row of the Stage 2 minimum-edit-set comparison table."""

    antigen: str
    sheep_to_baboon: EditNecessity
    pig_to_baboon: EditNecessity
    delta_note: str


class XenoEditPlan(_Frozen):
    """Stage 2 xenograft-back-to-baboon editing plan."""

    donor_species: Species
    recipient_species: Species
    rows: tuple[XenoEditRow, ...]
    guides: tuple[GuideRNA, ...]
    delta_vs_pig_protocol: str
    references: tuple[str, ...]


class VerificationResult(_Frozen):
    """Outcome of Tier 1 golden-standard verification on a single entry."""

    golden_id: str
    gene_symbol: str
    species: Species
    expected_protospacer: str
    found_rank: int | None  # None if not found in the scored output
    total_candidates: int
    passed: bool
