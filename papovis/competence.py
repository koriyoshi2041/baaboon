"""Stage 1 (donor side): competence-enhancing edits on baboon iPSCs.

Interspecies chimerism between primate and ungulate embryos is known to be
limited by at least two mechanistic classes (Zheng 2021, Lin 2024):
  - cell competition / apoptosis
  - cell-adhesion incompatibility

For each class this module emits a `CompetenceEditPlan` summarising:
  - which barrier gene is targeted
  - whether it is a knockout, knock-in, or assess-only candidate
  - the corresponding baboon guides (for KO targets)
  - a short rationale citing the literature

If Ensembl lacks a symbol annotation for a barrier gene in the donor species
(not uncommon for non-model primate genomes), the plan is still emitted
with zero guides and the rationale is extended to flag the missing annotation,
rather than failing the whole pipeline.
"""

from __future__ import annotations

import logging

from papovis.catalog import barrier_catalog, barrier_genes
from papovis.ensembl import EnsemblClient, EnsemblError
from papovis.grna import design_guides, top_n
from papovis.models import (
    CompetenceEditPlan,
    EditDirection,
    GuideRNA,
    Species,
)

logger = logging.getLogger(__name__)


def design_competence_edits(
    *,
    donor_species: Species = "Papio anubis",
    host_species: Species = "Ovis aries",
    client: EnsemblClient | None = None,
    guides_per_gene: int = 5,
    max_exons: int = 3,
) -> list[CompetenceEditPlan]:
    """Design iPSC competence-enhancing edits for every curated barrier gene."""
    client = client or EnsemblClient()
    plans: list[CompetenceEditPlan] = []
    catalog = barrier_catalog()
    classes = catalog.get("classes", {})

    for class_name, cls in classes.items():
        mechanism = cls.get("mechanism", "").strip()
        refs = tuple(cls.get("references", []))
        for gene in barrier_genes(class_name=class_name):
            direction = EditDirection(gene["direction"])
            symbol = gene["symbol"]
            rationale = (
                f"[{class_name}] "
                f"{gene.get('effect_summary') or gene.get('note', '')} "
                f"(mechanism: {mechanism})"
            ).strip()

            guides: tuple[GuideRNA, ...] = ()
            if direction is EditDirection.KO:
                try:
                    guides = tuple(
                        _design_guides_for_donor_gene(
                            client=client,
                            donor_species=donor_species,
                            gene_symbol=symbol,
                            guides_per_gene=guides_per_gene,
                            max_exons=max_exons,
                        )
                    )
                except EnsemblError as exc:
                    logger.warning(
                        "No Ensembl annotation for %s in %s; emitting plan "
                        "with zero guides. (%s)",
                        symbol,
                        donor_species,
                        exc,
                    )
                    rationale = (
                        f"{rationale} "
                        f"[NOTE: no Ensembl symbol annotation for {symbol} in "
                        f"{donor_species} at release 115; resolve via human "
                        f"ortholog lookup before wet-lab design.]"
                    )
            # OVEREXPRESS and ASSESS do not emit sgRNAs from this module.
            plans.append(
                CompetenceEditPlan(
                    donor_species=donor_species,
                    host_species=host_species,
                    barrier_class=class_name,
                    direction=direction,
                    gene_symbol=symbol,
                    guides=guides,
                    rationale=rationale,
                    references=refs,
                )
            )
    return plans


def _design_guides_for_donor_gene(
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
