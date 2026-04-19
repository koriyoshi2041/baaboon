"""MVP demo: baboon thymus in sheep, xenograft into baboon recipient.

Thymus is a compelling target because a donor-derived thymus engrafted in
the recipient educates recipient T cells in a donor-permissive way, which
can reduce rejection of subsequent solid organs from the same donor. The
niche gene is FOXN1, whose loss produces the athymic "nude" phenotype in
mouse/rat (Nehls 1994, Nature 372:103).
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
    path = Path("reports/thymus.md")

    t0 = time.time()
    client = EnsemblClient()

    print("Stage 1a — niche vacancy for thymus (FOXN1)")
    niche = design_niche_edits(
        organ="thymus",
        host_species="Ovis aries",
        client=client,
        guides_per_gene=10,
    )
    print(f"  {len(niche.guides)} guides against {', '.join(niche.target_genes)}")

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
