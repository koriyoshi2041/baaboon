"""Tier 3 biological-sanity check: phylogenetic identity ordering.

Human is closer to baboon than sheep is, so for any well-conserved gene
Ensembl homology should report ``human-baboon protein identity > sheep-baboon
protein identity``. If the ortholog resolver ever inverts this, something is
broken (wrong gene matched, assembly confusion, or simply unreliable data).

These tests are marked ``network``. Runs against rest.ensembl.org, response
cached under ~/.cache/papovis/.
"""

from __future__ import annotations

import pytest

from papovis.ensembl import EnsemblClient

# Well-conserved, protein-coding genes present in all three species. TP53
# and PDX1 are load-bearing in papovis itself; adding MSTN and FOXN1
# broadens the evidence base.
_CANONICAL_ORTHOLOGS = ["TP53", "PDX1", "MSTN", "FOXN1"]


@pytest.mark.network
@pytest.mark.parametrize("gene", _CANONICAL_ORTHOLOGS)
def test_primate_ortholog_closer_than_ungulate(gene: str) -> None:
    client = EnsemblClient()

    try:
        human_baboon = _best_protein_identity(
            client, source="Homo sapiens", target="Papio anubis", gene=gene
        )
        sheep_baboon = _best_protein_identity(
            client, source="Ovis aries", target="Papio anubis", gene=gene
        )
    except LookupError as exc:
        pytest.skip(str(exc))

    assert human_baboon > sheep_baboon, (
        f"Phylogeny violated for {gene}: "
        f"human-baboon identity {human_baboon:.1f}% is not greater than "
        f"sheep-baboon identity {sheep_baboon:.1f}%. "
        "Either the ortholog resolver is returning the wrong gene or the "
        "identity-parsing branch has regressed."
    )


def _best_protein_identity(
    client: EnsemblClient,
    *,
    source: str,
    target: str,
    gene: str,
) -> float:
    """Return max protein identity (source_species_perc_id vs target in %)."""
    homologies = client.fetch_orthologs(
        source_species=source,  # type: ignore[arg-type]
        gene_symbol=gene,
        target_species=target,  # type: ignore[arg-type]
    )
    if not homologies:
        raise LookupError(f"No {source}→{target} homology returned for {gene}")
    # Condensed format gives `type` and `target` stable ID but NOT percent
    # identity — fall back to the full format for that. Since our client
    # requests condensed by default, we treat the presence of any orthologue
    # as evidence of relatedness, and use sequence length as a proxy for
    # "close" matching when identity fields are unavailable.
    # Real identity requires format=full, which returns `perc_id` fields.
    full = client._get(  # type: ignore[attr-defined]
        f"/homology/symbol/{_species_to_name(source)}/{gene}",
        {
            "target_species": _species_to_name(target),
            "type": "orthologues",
            "format": "full",
        },
    )
    data = full.get("data", []) if isinstance(full, dict) else []
    if not data:
        raise LookupError(f"Empty data for {gene} {source}→{target}")
    identities: list[float] = []
    for entry in data[0].get("homologies", []):
        src = entry.get("source", {})
        perc_id = src.get("perc_id")
        if perc_id is not None:
            identities.append(float(perc_id))
    if not identities:
        raise LookupError(
            f"Homology payload for {gene} has no perc_id fields"
        )
    return max(identities)


def _species_to_name(latin: str) -> str:
    return {
        "Homo sapiens": "homo_sapiens",
        "Papio anubis": "papio_anubis",
        "Ovis aries": "ovis_aries",
    }[latin]
