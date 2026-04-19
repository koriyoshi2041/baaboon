"""Stage 2: minimum xeno-antigen edit set for sheep→baboon xenograft.

The headline insight of this module is the DELTA against the standard
pig→baboon (and pig→human) 10-gene-edit protocol exemplified by
United Therapeutics / eGenesis donor pigs (2024-2025). Two edits that are
obligate in pig donors become unnecessary or even counter-productive when
the donor is a sheep:

  * CMAH — pig lacks Neu5Gc preservation relative to human, but baboon
    retains CMAH; sheep also retains CMAH. Estrada 2015 showed CMAH KO in
    pig donors INCREASES baboon serum binding 3-fold.

  * (Candidate) B4GALNT2 — baboon has partial Sda expression; wet-lab
    confirmation needed on sheep. Marked RECOMMENDED rather than REQUIRED.

This module compiles the side-by-side edit table and, where a KO is still
required, designs sgRNAs against the sheep genome. For knock-in targets
(complement regulators CD46/CD55/CD59, thrombomodulin, EPCR) the module
annotates the donor cDNA source as the baboon (not human) ortholog.
"""

from __future__ import annotations

import logging

from papovis.catalog import xeno_antigens
from papovis.ensembl import EnsemblClient, EnsemblError
from papovis.grna import design_guides, top_n
from papovis.models import (
    EditNecessity,
    GuideRNA,
    Species,
    XenoEditPlan,
    XenoEditRow,
)

logger = logging.getLogger(__name__)


def design_xeno_edits(
    *,
    donor_species: Species = "Ovis aries",
    recipient_species: Species = "Papio anubis",
    client: EnsemblClient | None = None,
    guides_per_gene: int = 8,
    max_exons: int = 3,
) -> XenoEditPlan:
    """Return the Stage 2 minimum-edit-set plan for sheep→baboon xenograft."""
    client = client or EnsemblClient()
    catalog = xeno_antigens()
    rows: list[XenoEditRow] = []
    guides: list[GuideRNA] = []
    references: list[str] = []

    # Glycan antigens — the only group where the sheep/pig delta emerges.
    for ag in catalog["glycan_antigens"]:
        sheep_flag = EditNecessity(ag["sheep_to_baboon_edit"])
        pig_flag = EditNecessity(ag["pig_to_baboon_edit"])
        rows.append(
            XenoEditRow(
                antigen=ag["symbol"],
                sheep_to_baboon=sheep_flag,
                pig_to_baboon=pig_flag,
                delta_note=str(ag.get("delta_note", "")),
            )
        )
        references.extend(ag.get("references", []))
        if sheep_flag in (EditNecessity.REQUIRED, EditNecessity.RECOMMENDED):
            try:
                guides.extend(
                    _design_sheep_guides(
                        client=client,
                        donor_species=donor_species,
                        gene_symbol=ag["symbol"],
                        guides_per_gene=guides_per_gene,
                        max_exons=max_exons,
                    )
                )
            except EnsemblError as exc:
                logger.warning(
                    "Ensembl symbol %s not directly resolvable in %s; "
                    "wet-lab curators should cross-reference via ortholog "
                    "mapping (NCBI Gene, UniProt, or BioMart). (%s)",
                    ag["symbol"],
                    donor_species,
                    exc,
                )

    # Knock-in (complement + coagulation) — annotated as KI rows with zero
    # sgRNAs emitted here; the cDNA source selection (Papio anubis) is the
    # relevant decision.
    for ki in [*catalog["complement_regulators"], *catalog["coagulation_thrombosis"]]:
        rows.append(
            XenoEditRow(
                antigen=ki["symbol"],
                sheep_to_baboon=EditNecessity(ki["sheep_to_baboon_edit"]),
                pig_to_baboon=EditNecessity.REQUIRED,
                delta_note=f"Knock-in. Donor cDNA source: {ki['donor_cdna_source']}",
            )
        )

    delta = _compute_delta(rows)
    return XenoEditPlan(
        donor_species=donor_species,
        recipient_species=recipient_species,
        rows=tuple(rows),
        guides=tuple(guides),
        delta_vs_pig_protocol=delta,
        references=tuple(dict.fromkeys(references)),
    )


def _design_sheep_guides(
    *,
    client: EnsemblClient,
    donor_species: Species,
    gene_symbol: str,
    guides_per_gene: int,
    max_exons: int,
) -> list[GuideRNA]:
    exons = client.fetch_exon_sequences(
        species=donor_species,
        gene_symbol=gene_symbol,
        max_exons=max_exons,
    )
    per_gene: list[GuideRNA] = []
    for exon in exons:
        per_gene.extend(
            design_guides(
                sequence=exon.sequence,
                species=donor_species,
                gene_symbol=gene_symbol,
                chromosome=exon.chromosome,
                region_start=exon.start,
                region_strand=exon.strand,
            )
        )
    per_gene.sort(key=lambda g: (-g.on_target_score, g.start))
    return top_n(per_gene, guides_per_gene)


def _compute_delta(rows: list[XenoEditRow]) -> str:
    pig_required = sum(1 for r in rows if r.pig_to_baboon is EditNecessity.REQUIRED)
    sheep_required = sum(1 for r in rows if r.sheep_to_baboon is EditNecessity.REQUIRED)
    sheep_counter = sum(
        1 for r in rows if r.sheep_to_baboon is EditNecessity.NOT_REQUIRED
    )
    return (
        f"Pig→baboon REQUIRED edits: {pig_required}. "
        f"Sheep→baboon REQUIRED edits: {sheep_required}. "
        f"Edits DROPPED on sheep axis: {sheep_counter}."
    )
