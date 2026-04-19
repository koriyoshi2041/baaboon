"""Typer-based CLI. Example:

    papovis plan --organ pancreas --output reports/pancreas.md
    papovis organs
    papovis verify --gene PDX1 --species "Ovis aries"
"""

from __future__ import annotations

from pathlib import Path

import typer
from rich.console import Console

from papovis.catalog import available_organs
from papovis.competence import design_competence_edits
from papovis.design import design_guides_for_gene
from papovis.ensembl import EnsemblClient
from papovis.golden import summarise, verify_guides_against_golden
from papovis.niche import design_niche_edits
from papovis.report import render_markdown
from papovis.xeno import design_xeno_edits

app = typer.Typer(
    add_completion=False,
    help="papovis — two-stage CRISPR design for baboon×sheep chimera + xenograft.",
)
console = Console()


@app.command()
def organs() -> None:
    """List organs for which Stage-1 niche genes are curated."""
    for organ in available_organs():
        console.print(f"- {organ}")


@app.command()
def plan(
    organ: str = typer.Option(..., help="Target organ, e.g. pancreas / kidney / thymus."),
    output: Path = typer.Option(Path("reports/plan.md"), help="Markdown output path."),
    guides_per_gene: int = typer.Option(10, min=1, max=50),
) -> None:
    """Build a full two-stage plan and write a markdown report."""
    client = EnsemblClient()
    console.print(f"[bold]Stage 1a[/bold] — niche edits in sheep for {organ}")
    niche = design_niche_edits(
        organ=organ, client=client, guides_per_gene=guides_per_gene
    )
    console.print(f"[bold]Stage 1b[/bold] — competence edits on baboon iPSCs")
    competence = design_competence_edits(client=client)
    console.print("[bold]Stage 2[/bold] — sheep→baboon xeno-antigen minimum set")
    xeno = design_xeno_edits(client=client)

    markdown = render_markdown(niche=niche, competence=competence, xeno=xeno)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(markdown, encoding="utf-8")
    console.print(f"[green]Report written to[/green] {output}")


@app.command()
def verify(
    gene: str = typer.Option(..., help="Gene symbol to verify, e.g. PDX1."),
    species: str = typer.Option("Ovis aries", help="Latin species name."),
    guides_per_gene: int = typer.Option(20, min=1, max=100),
    top_n: int = typer.Option(20, min=1, max=100, help="Pass threshold: rank ≤ top_n."),
) -> None:
    """Run Tier-1 golden-standard verification for any gene.

    Fetches the gene from Ensembl, designs sgRNAs, and compares to every
    entry in data/golden_standards.yaml that matches (gene, species).
    """
    client = EnsemblClient()
    guides = design_guides_for_gene(
        gene_symbol=gene,
        species=species,  # type: ignore[arg-type]
        client=client,
        guides_per_gene=guides_per_gene,
    )
    results = verify_guides_against_golden(
        produced_guides=guides,
        gene_symbol=gene,
        species=species,  # type: ignore[arg-type]
        top_n_for_pass=top_n,
    )
    console.print(summarise(results))


if __name__ == "__main__":
    app()
