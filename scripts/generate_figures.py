"""Run the full papovis pipeline and produce publication-quality figures.

Outputs:
    figures/fig1_edit_delta.{pdf,png}   — Stage-2 edit burden (project thesis)
    figures/fig2_verification.{pdf,png} — Tier-1 MSTN + PDX1 recapitulation
    figures/fig3_cross_organ.{pdf,png}  — Niche sgRNA landscape across organs

All plots are generated from REAL Ensembl REST output (cached under
~/.cache/papovis). Re-running the script costs one HTTP round-trip per
gene the first time; subsequent runs reuse the cache.
"""

from __future__ import annotations

import logging
import time
from pathlib import Path

from papovis.catalog import available_organs, niche_genes_for_organ
from papovis.design import design_guides_for_gene
from papovis.ensembl import EnsemblClient, EnsemblError
from papovis.figures import (
    figure_cross_organ_landscape,
    figure_edit_delta,
    figure_verification,
    save,
)
from papovis.xeno import design_xeno_edits

FIG_DIR = Path("figures")
MSTN_PROTOSPACER = "GGCTGTGTAATGCATGCTTG"
PDX1_SINGLE_PROTOSPACER = "GGGCCCCGCTGGAACGCGCA"


def main() -> None:
    logging.basicConfig(level=logging.WARNING, format="  [warn] %(message)s")
    FIG_DIR.mkdir(exist_ok=True)
    t0 = time.time()
    client = EnsemblClient()

    print("=== Figure 1: Stage-2 edit delta ===")
    xeno = design_xeno_edits(client=client, guides_per_gene=5)
    fig1 = figure_edit_delta(xeno)
    save(fig1, FIG_DIR / "fig1_edit_delta")
    print(f"  [{time.time()-t0:.1f}s] wrote fig1_edit_delta.{{pdf,png}}")
    print(f"  delta: {xeno.delta_vs_pig_protocol}")

    print("=== Figure 2: Tier-1 verification ===")
    mstn_guides = design_guides_for_gene(
        gene_symbol="MSTN", species="Ovis aries",
        client=client, guides_per_gene=500,
    )
    pdx1_guides = design_guides_for_gene(
        gene_symbol="PDX1", species="Ovis aries",
        client=client, guides_per_gene=500,
    )
    fig2 = figure_verification(
        mstn_guides=mstn_guides,
        pdx1_guides=pdx1_guides,
        mstn_published=MSTN_PROTOSPACER,
        pdx1_published=PDX1_SINGLE_PROTOSPACER,
    )
    save(fig2, FIG_DIR / "fig2_verification")
    print(f"  [{time.time()-t0:.1f}s] wrote fig2_verification.{{pdf,png}}")
    print(f"  MSTN {len(mstn_guides)} guides, PDX1 {len(pdx1_guides)} guides")

    print("=== Figure 3: Cross-organ niche landscape ===")
    organ_to_gene = {}
    for organ in available_organs():
        inner = {}
        for gene_entry in niche_genes_for_organ(organ):
            sym = gene_entry["symbol"]
            try:
                inner[sym] = design_guides_for_gene(
                    gene_symbol=sym, species="Ovis aries",
                    client=client, guides_per_gene=2000,
                )
            except EnsemblError as exc:
                logging.warning("skip %s/%s: %s", organ, sym, exc)
        organ_to_gene[organ] = inner
    fig3 = figure_cross_organ_landscape(organ_to_gene)
    save(fig3, FIG_DIR / "fig3_cross_organ")
    print(f"  [{time.time()-t0:.1f}s] wrote fig3_cross_organ.{{pdf,png}}")

    print(f"\nAll figures in {FIG_DIR.resolve()} ({time.time()-t0:.1f}s total)")


if __name__ == "__main__":
    main()
