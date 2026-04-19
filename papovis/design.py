"""Generic, catalog-free gRNA design for arbitrary gene symbols.

Used by the `verify` CLI path and by tests that need to check recovery of
literature sgRNAs on loci that aren't niche or xeno genes (e.g. MSTN).
"""

from __future__ import annotations

from papovis.ensembl import EnsemblClient
from papovis.grna import design_guides, top_n
from papovis.models import GuideRNA, Species


def design_guides_for_gene(
    *,
    gene_symbol: str,
    species: Species,
    client: EnsemblClient | None = None,
    guides_per_gene: int = 20,
    max_exons: int = 3,
    pam: str = "NGG",
    guide_length: int = 20,
) -> list[GuideRNA]:
    """Fetch the first N exons of a gene and return ranked sgRNAs."""
    client = client or EnsemblClient()
    exons = client.fetch_exon_sequences(
        species=species,
        gene_symbol=gene_symbol,
        max_exons=max_exons,
    )
    all_guides: list[GuideRNA] = []
    for exon in exons:
        all_guides.extend(
            design_guides(
                sequence=exon.sequence,
                species=species,
                gene_symbol=gene_symbol,
                chromosome=exon.chromosome,
                region_start=exon.start,
                region_strand=exon.strand,
                pam=pam,
                guide_length=guide_length,
            )
        )
    all_guides.sort(key=lambda g: (-g.on_target_score, g.start))
    return top_n(all_guides, guides_per_gene)
