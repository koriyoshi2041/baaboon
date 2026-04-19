"""Tier 3 biological-sanity check on the Stage 2 plan.

These assertions encode the *project thesis*: that choosing sheep rather than
pig as the donor to a baboon recipient reduces the obligate edit count by at
least one (CMAH). If the xeno_antigens.yaml is edited in a way that erases
this differential, CI must fail.
"""

from __future__ import annotations

from papovis.catalog import xeno_antigens
from papovis.models import EditNecessity


def test_cmah_delta_is_preserved() -> None:
    glycans = xeno_antigens()["glycan_antigens"]
    cmah = next(g for g in glycans if g["symbol"] == "CMAH")
    assert (
        cmah["sheep_to_baboon_edit"] == EditNecessity.NOT_REQUIRED.value
    ), "CMAH must remain NOT_REQUIRED for sheep→baboon (Estrada 2015)."


def test_ggta1_remains_required() -> None:
    glycans = xeno_antigens()["glycan_antigens"]
    ggta1 = next(g for g in glycans if g["symbol"] == "GGTA1")
    assert (
        ggta1["sheep_to_baboon_edit"] == EditNecessity.REQUIRED.value
    ), "GGTA1 must remain REQUIRED on sheep→baboon axis (Galili 1988)."


def test_aggregate_delta_is_positive() -> None:
    """The number of edits dropped on the sheep axis must be >= 1 for the
    project's central hypothesis to hold."""
    glycans = xeno_antigens()["glycan_antigens"]
    dropped = sum(
        1
        for g in glycans
        if g["sheep_to_baboon_edit"] == EditNecessity.NOT_REQUIRED.value
        and g["pig_to_baboon_edit"]
        in (
            EditNecessity.REQUIRED.value,
            EditNecessity.COUNTERPRODUCTIVE.value,
        )
    )
    assert dropped >= 1, (
        "Sheep→baboon hypothesis requires at least one edit dropped vs pig "
        "protocol; current YAML drops zero — either the curation is wrong or "
        "the project thesis needs revisiting."
    )
