"""Tier 1 verification: confirm that published sgRNAs from the literature
appear in the top-N output of our pipeline. Unpopulated entries (pending
manual curation) are returned as `VerificationResult(passed=False,
found_rank=None)` so the verification harness can distinguish a "not yet
curated" state from a real failure.
"""

from __future__ import annotations

from papovis.catalog import golden_entries
from papovis.models import GuideRNA, Species, VerificationResult


def verify_guides_against_golden(
    *,
    produced_guides: list[GuideRNA],
    gene_symbol: str,
    species: Species,
    top_n_for_pass: int = 10,
) -> list[VerificationResult]:
    """Compare produced guides to every golden-standard entry for (gene, species)."""
    results: list[VerificationResult] = []
    relevant_entries = [
        e
        for e in golden_entries()
        if e.get("gene_symbol") == gene_symbol and e.get("species") == species
    ]

    for entry in relevant_entries:
        expected = entry.get("protospacer")
        total = len(produced_guides)
        allow_snp_miss = bool(entry.get("allow_snp_miss"))
        if not expected:
            results.append(
                VerificationResult(
                    golden_id=entry["id"],
                    gene_symbol=gene_symbol,
                    species=species,
                    expected_protospacer="",
                    found_rank=None,
                    total_candidates=total,
                    passed=False,
                )
            )
            continue
        expected_upper = expected.upper()
        rank: int | None = None
        for idx, guide in enumerate(produced_guides, start=1):
            if guide.protospacer.upper() == expected_upper:
                rank = idx
                break
        # SNP-tolerant entries pass by default when not found (because the
        # current reference genome is known to differ from the paper's).
        if rank is None and allow_snp_miss:
            passed = True
        else:
            passed = rank is not None and rank <= top_n_for_pass
        results.append(
            VerificationResult(
                golden_id=entry["id"],
                gene_symbol=gene_symbol,
                species=species,
                expected_protospacer=expected_upper,
                found_rank=rank,
                total_candidates=total,
                passed=passed,
            )
        )

    return results


def summarise(results: list[VerificationResult]) -> str:
    if not results:
        return "No golden-standard entries for this target."
    pending = [r for r in results if not r.expected_protospacer]
    scored = [r for r in results if r.expected_protospacer]
    lines = []
    if pending:
        lines.append(
            f"⏳ {len(pending)} golden entry/entries pending manual curation "
            "(protospacer sequence empty in golden_standards.yaml)."
        )
    for r in scored:
        state = "✅" if r.passed else "❌"
        rank_str = str(r.found_rank) if r.found_rank is not None else "not found"
        lines.append(
            f"{state} {r.golden_id}: rank={rank_str} "
            f"(of {r.total_candidates}) for {r.gene_symbol} in {r.species}."
        )
    return "\n".join(lines)
