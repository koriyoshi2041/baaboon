"""Tier-1 live verification against Crispo et al. 2015 PLOS ONE MSTN sheep.

Marked ``network`` because it hits `rest.ensembl.org`. Run with:

    pytest -m network tests/test_crispo_mstn.py -v

Expected outcome: the published protospacer ``GGCTGTGTAATGCATGCTTG`` with
PAM ``TGG`` appears at the forward strand of sheep MSTN exon 1 and is
ranked ≤ 20 by our scorer.
"""

from __future__ import annotations

import pytest

from papovis.design import design_guides_for_gene
from papovis.golden import verify_guides_against_golden


@pytest.mark.golden
@pytest.mark.network
def test_crispo_mstn_protospacer_is_recovered() -> None:
    guides = design_guides_for_gene(
        gene_symbol="MSTN",
        species="Ovis aries",
        guides_per_gene=50,
    )
    assert guides, "No guides returned for sheep MSTN"

    results = verify_guides_against_golden(
        produced_guides=guides,
        gene_symbol="MSTN",
        species="Ovis aries",
        top_n_for_pass=20,
    )
    # Only evaluate entries that have a populated protospacer.
    populated = [r for r in results if r.expected_protospacer]
    assert populated, (
        "Expected at least one curated golden-standard entry for MSTN sheep."
    )
    failed = [r for r in populated if not r.passed]
    assert not failed, (
        "Pipeline did not recover Crispo 2015 MSTN sgRNA in top-20: "
        + ", ".join(
            f"{r.golden_id} (found_rank={r.found_rank})" for r in failed
        )
    )
