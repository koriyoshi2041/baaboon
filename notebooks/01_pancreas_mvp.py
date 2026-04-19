"""MVP demo: design the full two-stage plan for a baboon-derived pancreas
grown in a sheep, then xeno-transplanted into a baboon recipient.

Runs against live Ensembl REST (release 115) and writes a markdown report.

Usage:
    python notebooks/01_pancreas_mvp.py

Requires network access to rest.ensembl.org. Responses are cached under
`~/.cache/papovis/` so subsequent runs are instantaneous.
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
    logging.basicConfig(
        level=logging.WARNING,
        format="  [warn] %(message)s",
    )
    report_dir = Path("reports")
    report_dir.mkdir(exist_ok=True)
    report_path = report_dir / "pancreas.md"

    t0 = time.time()
    client = EnsemblClient()

    print("Stage 1a — host-embryo niche vacancy in sheep")
    niche = design_niche_edits(
        organ="pancreas",
        host_species="Ovis aries",
        client=client,
        guides_per_gene=10,
    )
    print(f"  {len(niche.guides)} guides across {len(niche.target_genes)} gene(s)")

    print("Stage 1b — chimera-barrier edits on baboon iPSCs")
    competence = design_competence_edits(client=client, guides_per_gene=3)
    with_guides = sum(1 for p in competence if p.guides)
    print(f"  {len(competence)} plans ({with_guides} with KO guides emitted)")

    print("Stage 2 — sheep→baboon xeno-antigen minimum edit set")
    xeno = design_xeno_edits(client=client, guides_per_gene=5)
    print(f"  {len(xeno.guides)} sheep-side KO guides")
    print(f"  DELTA: {xeno.delta_vs_pig_protocol}")

    markdown = render_markdown(niche=niche, competence=competence, xeno=xeno)
    report_path.write_text(markdown, encoding="utf-8")
    print(f"\nReport written to {report_path} ({len(markdown)} chars)")
    print(f"Elapsed: {time.time() - t0:.1f}s")

    _print_headline(niche=niche, xeno=xeno)


def _print_headline(*, niche, xeno) -> None:  # type: ignore[no-untyped-def]
    """Emit the single-line headline that verifies the project thesis."""
    top_guide = niche.guides[0] if niche.guides else None
    print()
    print("=" * 72)
    print("HEADLINE")
    print("=" * 72)
    if top_guide:
        print(
            f"Top PDX1 sheep sgRNA: {top_guide.protospacer}{top_guide.pam} "
            f"(chr{top_guide.chromosome}:{top_guide.start}, strand {top_guide.strand}, "
            f"score {top_guide.on_target_score:.3f})"
        )
    print(f"Stage 2 edit delta: {xeno.delta_vs_pig_protocol}")
    print("=" * 72)


if __name__ == "__main__":
    main()
