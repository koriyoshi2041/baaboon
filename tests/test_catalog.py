"""Validate that the curated YAML catalogs have the structure downstream
modules assume. These tests run offline and guard against accidental edits.
"""

from __future__ import annotations

import pytest

from papovis.catalog import (
    available_organs,
    barrier_genes,
    golden_entries,
    niche_genes_for_organ,
    xeno_antigens,
)
from papovis.models import EditDirection, EditNecessity


class TestNicheCatalog:
    def test_pancreas_present(self) -> None:
        assert "pancreas" in available_organs()

    def test_pancreas_uses_pdx1(self) -> None:
        genes = niche_genes_for_organ("pancreas")
        symbols = [g["symbol"] for g in genes]
        assert "PDX1" in symbols

    def test_kidney_uses_sall1(self) -> None:
        genes = niche_genes_for_organ("kidney")
        symbols = [g["symbol"] for g in genes]
        assert "SALL1" in symbols

    def test_all_niche_genes_have_references(self) -> None:
        for organ in available_organs():
            for gene in niche_genes_for_organ(organ):
                assert gene.get("references"), (
                    f"{gene['symbol']} (organ={organ}) is missing references; "
                    "every curated niche gene must cite its primary literature."
                )


class TestBarrierCatalog:
    def test_two_barrier_classes(self) -> None:
        classes = {g["barrier_class"] for g in barrier_genes()}
        assert "apoptosis_cell_competition" in classes
        assert "cell_adhesion" in classes

    def test_tp53_is_ko(self) -> None:
        tp53 = [g for g in barrier_genes() if g["symbol"] == "TP53"]
        assert tp53 and tp53[0]["direction"] == EditDirection.KO.value

    def test_bcl2_is_overexpress(self) -> None:
        bcl2 = [g for g in barrier_genes() if g["symbol"] == "BCL2"]
        assert bcl2 and bcl2[0]["direction"] == EditDirection.OVEREXPRESS.value

    def test_cdh1_is_assess_only(self) -> None:
        cdh1 = [g for g in barrier_genes() if g["symbol"] == "CDH1"]
        assert cdh1 and cdh1[0]["direction"] == EditDirection.ASSESS.value


class TestXenoCatalog:
    def test_cmah_is_not_required_for_sheep_baboon(self) -> None:
        """Tier 3 sanity check: reflects Estrada 2015 finding."""
        glycans = xeno_antigens()["glycan_antigens"]
        cmah = next(g for g in glycans if g["symbol"] == "CMAH")
        assert (
            cmah["sheep_to_baboon_edit"] == EditNecessity.NOT_REQUIRED.value
        )
        assert (
            cmah["pig_to_baboon_edit"] == EditNecessity.COUNTERPRODUCTIVE.value
        )

    def test_ggta1_still_required_on_sheep_axis(self) -> None:
        glycans = xeno_antigens()["glycan_antigens"]
        ggta1 = next(g for g in glycans if g["symbol"] == "GGTA1")
        assert ggta1["sheep_to_baboon_edit"] == EditNecessity.REQUIRED.value

    def test_complement_regulators_use_baboon_cdna(self) -> None:
        for ki in xeno_antigens()["complement_regulators"]:
            assert "Papio" in ki["donor_cdna_source"], (
                f"{ki['symbol']} should use baboon cDNA for knock-in "
                "(the recipient species), not human."
            )


class TestGoldenStandards:
    def test_pdx1_entries_present(self) -> None:
        pdx1 = [e for e in golden_entries() if e["gene_symbol"] == "PDX1"]
        assert len(pdx1) >= 1, "Expected at least one PDX1 golden entry"

    def test_every_entry_has_doi(self) -> None:
        for entry in golden_entries():
            paper = entry.get("paper") or {}
            assert paper.get("doi"), (
                f"Golden entry {entry['id']} is missing a DOI; every "
                "reference protospacer must trace to a citable source."
            )

    def test_pending_curation_entries_have_null_sequence(self) -> None:
        for entry in golden_entries():
            if entry.get("status") == "pending_curation":
                assert entry.get("protospacer") in (None, ""), (
                    f"Entry {entry['id']} is marked pending_curation but "
                    "has a non-empty protospacer — update status to 'curated'."
                )


if __name__ == "__main__":  # pragma: no cover
    pytest.main([__file__, "-v"])
