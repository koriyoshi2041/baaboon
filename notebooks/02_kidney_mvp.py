"""MVP demo: design a full two-stage plan for growing a baboon-derived
kidney inside a sheep, then transplanting it back into a baboon recipient.

The niche-vacancy step targets SALL1 and SIX1 (rather than PDX1), following
the Wang et al. 2023 Cell Stem Cell pig-humanised-mesonephros precedent.
Runs against live Ensembl REST and writes reports/kidney.md.
"""

from __future__ import annotations

import logging
import time
from pathlib import Path

from papovis.competence import design_competence_edits
from papovis.ensembl import EnsemblClient
from papovis.niche import design_niche_edits
from papovis.report import render_markdown
from papovis.xeno import design_xeno_edits


def main() -> None:
    logging.basicConfig(level=logging.WARNING, format="  [warn] %(message)s")
    Path("reports").mkdir(exist_ok=True)
    path = Path("reports/kidney.md")

    t0 = time.time()
    client = EnsemblClient()

    print("Stage 1a — niche vacancy for kidney (SALL1 + SIX1)")
    niche = design_niche_edits(
        organ="kidney",
        host_species="Ovis aries",
        client=client,
        guides_per_gene=10,
    )
    print(f"  {len(niche.guides)} guides across {len(niche.target_genes)} gene(s)")

    print("Stage 1b — baboon-iPSC chimera-barrier edits")
    competence = design_competence_edits(client=client, guides_per_gene=3)

    print("Stage 2 — sheep→baboon xeno antigens")
    xeno = design_xeno_edits(client=client, guides_per_gene=5)
    print(f"  DELTA: {xeno.delta_vs_pig_protocol}")

    path.write_text(
        render_markdown(niche=niche, competence=competence, xeno=xeno),
        encoding="utf-8",
    )
    print(f"\nReport: {path}")
    print(f"Elapsed: {time.time() - t0:.1f}s")


if __name__ == "__main__":
    main()
