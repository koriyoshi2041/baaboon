"""Stage 1: niche-vacancy CRISPR design for the sheep host embryo.

Given a target organ (pancreas, kidney, thymus, ...), the niche module:

1. reads the curated `niche_genes.yaml` to learn the master developmental
   regulator(s) whose knockout creates the required developmental niche;
2. fetches the first few coding exons of each gene from Ensembl in the
   sheep (Ovis aries) genome;
3. scans each exon for SpCas9 guides and scores them;
4. emits a `NicheEditPlan` combining all genes and their top guides.

The design deliberately targets early coding exons rather than 5'UTR or
promoter because a frameshift in an early CDS reliably ablates protein
production, matching the experimental strategy of Vilarino 2017 (PDX1).
"""

from __future__ import annotations

from typing import Any

from papovis.catalog import niche_genes_for_organ
from papovis.ensembl import EnsemblClient
from papovis.grna import design_guides, top_n
from papovis.models import GuideRNA, NicheEditPlan, Species


def design_niche_edits(
    *,
    organ: str,
    host_species: Species = "Ovis aries",
    client: EnsemblClient | None = None,
    guides_per_gene: int = 10,
    max_exons: int = 3,
) -> NicheEditPlan:
    """Build the Stage 1 niche-vacancy CRISPR plan for one organ in the host."""
    genes = niche_genes_for_organ(organ)
    if not genes:
        raise ValueError(f"No niche genes curated for organ {organ!r}")

    client = client or EnsemblClient()
    all_guides: list[GuideRNA] = []
    references: list[str] = []
    target_symbols: list[str] = []

    for gene_entry in genes:
        symbol = gene_entry["symbol"]
        target_symbols.append(symbol)
        references.extend(gene_entry.get("references", []))
        exons = client.fetch_exon_sequences(
            species=host_species,
            gene_symbol=symbol,
            max_exons=max_exons,
        )
        per_gene: list[GuideRNA] = []
        for exon in exons:
            exon_guides = design_guides(
                sequence=exon.sequence,
                species=host_species,
                gene_symbol=symbol,
                chromosome=exon.chromosome,
                region_start=exon.start,
                region_strand=exon.strand,
            )
            per_gene.extend(exon_guides)
        per_gene.sort(key=lambda g: (-g.on_target_score, g.start))
        all_guides.extend(top_n(per_gene, guides_per_gene))

    return NicheEditPlan(
        organ=organ,
        host_species=host_species,
        target_genes=tuple(target_symbols),
        guides=tuple(all_guides),
        references=tuple(_dedupe(references)),
    )


def _dedupe(items: list[Any]) -> list[Any]:
    seen: set[Any] = set()
    out: list[Any] = []
    for item in items:
        if item not in seen:
            seen.add(item)
            out.append(item)
    return out
