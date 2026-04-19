"""Emit a human-readable markdown report combining Stage 1 + Stage 2 plans."""

from __future__ import annotations

from io import StringIO

from papovis.models import (
    CompetenceEditPlan,
    EditNecessity,
    GuideRNA,
    NicheEditPlan,
    XenoEditPlan,
)


def render_markdown(
    *,
    niche: NicheEditPlan,
    competence: list[CompetenceEditPlan],
    xeno: XenoEditPlan,
) -> str:
    buf = StringIO()
    _header(buf, niche=niche, xeno=xeno)
    _stage1_niche(buf, niche)
    _stage1_competence(buf, competence)
    _stage2_xeno(buf, xeno)
    _references(buf, niche=niche, competence=competence, xeno=xeno)
    return buf.getvalue()


# ---------------------------------------------------------------- sections


def _header(buf: StringIO, *, niche: NicheEditPlan, xeno: XenoEditPlan) -> None:
    buf.write(f"# papovis report — {niche.organ}\n\n")
    buf.write(
        f"**Host embryo:** {niche.host_species}  \n"
        f"**Donor iPSCs:** {xeno.recipient_species}  \n"
        f"**Xenograft recipient:** {xeno.recipient_species}\n\n"
    )
    buf.write(
        "This report is a computational design artefact. It is not a wet-lab "
        "protocol and it does not supersede peer-reviewed publications.\n\n"
    )


def _stage1_niche(buf: StringIO, plan: NicheEditPlan) -> None:
    buf.write("## Stage 1a — host-embryo niche vacancy\n\n")
    buf.write(
        f"Knocking out {', '.join(plan.target_genes)} in a "
        f"{plan.host_species} zygote creates a developmental niche for the "
        f"{plan.organ}. Top guides:\n\n"
    )
    _guide_table(buf, plan.guides)


def _stage1_competence(buf: StringIO, plans: list[CompetenceEditPlan]) -> None:
    buf.write("\n## Stage 1b — donor-iPSC chimera-barrier edits\n\n")
    if not plans:
        buf.write("_No competence edits configured._\n")
        return
    for plan in plans:
        buf.write(
            f"### {plan.gene_symbol}  \n"
            f"Direction: **{plan.direction.value}**  \n"
            f"Class: {plan.barrier_class}  \n"
            f"Rationale: {plan.rationale}\n\n"
        )
        if plan.guides:
            _guide_table(buf, plan.guides)
        else:
            buf.write("_No sgRNAs emitted for this direction._\n\n")


def _stage2_xeno(buf: StringIO, plan: XenoEditPlan) -> None:
    buf.write("\n## Stage 2 — sheep→baboon xeno-antigen minimum edit set\n\n")
    buf.write(f"**Delta vs pig protocol:** {plan.delta_vs_pig_protocol}\n\n")
    buf.write("| Antigen | Sheep→Baboon | Pig→Baboon | Note |\n")
    buf.write("| --- | --- | --- | --- |\n")
    for row in plan.rows:
        buf.write(
            f"| {row.antigen} | {_fmt_necessity(row.sheep_to_baboon)} | "
            f"{_fmt_necessity(row.pig_to_baboon)} | {row.delta_note} |\n"
        )
    buf.write("\n")
    if plan.guides:
        buf.write("### Sheep-side knock-out sgRNAs\n\n")
        _guide_table(buf, plan.guides)


def _references(
    buf: StringIO,
    *,
    niche: NicheEditPlan,
    competence: list[CompetenceEditPlan],
    xeno: XenoEditPlan,
) -> None:
    buf.write("\n## References\n\n")
    seen: set[str] = set()
    for ref in list(niche.references) + list(xeno.references):
        if ref not in seen:
            seen.add(ref)
            buf.write(f"- {ref}\n")
    for plan in competence:
        for ref in plan.references:
            if ref not in seen:
                seen.add(ref)
                buf.write(f"- {ref}\n")


# ----------------------------------------------------------------- helpers


def _guide_table(buf: StringIO, guides: tuple[GuideRNA, ...]) -> None:
    if not guides:
        buf.write("_No guides emitted._\n")
        return
    buf.write(
        "| Gene | Protospacer (5'→3') | PAM | Chrom | Pos | Strand | GC | Score |\n"
        "| --- | --- | --- | --- | --- | --- | --- | --- |\n"
    )
    for g in guides:
        buf.write(
            f"| {g.gene_symbol} | `{g.protospacer}` | {g.pam} | "
            f"{g.chromosome} | {g.start} | {g.strand} | {g.gc_content:.2f} | "
            f"{g.on_target_score:.3f} |\n"
        )


def _fmt_necessity(value: EditNecessity) -> str:
    markers = {
        EditNecessity.REQUIRED: "**REQUIRED**",
        EditNecessity.RECOMMENDED: "recommended",
        EditNecessity.NOT_REQUIRED: "~~dropped~~",
        EditNecessity.COUNTERPRODUCTIVE: "⚠ counterproductive",
        EditNecessity.ASSESS: "assess",
    }
    return markers.get(value, value.value)
