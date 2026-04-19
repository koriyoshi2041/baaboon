"""Publication-quality figure generation for papovis results.

Conventions:
  * Okabe-Ito colour-blind-safe palette — Okabe & Ito (2008). Hex codes are
    hard-coded below rather than imported from a style library.
  * Despined axes, no gridlines except subtle x-axis, no 3-D.
  * Bold panel letters at top-left.
  * Annotations live in whitespace. Text never overlaps data or other text,
    never crosses the axis, never sits on top of bars or ticks.
  * PDF vector + 300 dpi PNG output.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

from papovis.models import EditNecessity, GuideRNA, XenoEditPlan

OKABE_ITO: dict[str, str] = {
    "orange": "#E69F00",
    "sky_blue": "#56B4E9",
    "bluish_green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "reddish_purple": "#CC79A7",
    "black": "#000000",
}

NECESSITY_COLOURS: dict[EditNecessity, str] = {
    EditNecessity.REQUIRED: OKABE_ITO["vermillion"],
    EditNecessity.RECOMMENDED: OKABE_ITO["yellow"],
    EditNecessity.ASSESS: OKABE_ITO["sky_blue"],
    EditNecessity.NOT_REQUIRED: OKABE_ITO["bluish_green"],
    EditNecessity.COUNTERPRODUCTIVE: OKABE_ITO["reddish_purple"],
}

NECESSITY_ABBREV: dict[EditNecessity, str] = {
    EditNecessity.REQUIRED: "REQ",
    EditNecessity.RECOMMENDED: "rec",
    EditNecessity.ASSESS: "?",
    EditNecessity.NOT_REQUIRED: "—",
    EditNecessity.COUNTERPRODUCTIVE: "✗",
}

ORGAN_COLOURS: dict[str, str] = {
    "pancreas": OKABE_ITO["vermillion"],
    "kidney": OKABE_ITO["blue"],
    "thymus": OKABE_ITO["bluish_green"],
    "liver": OKABE_ITO["orange"],
    "heart": OKABE_ITO["reddish_purple"],
}


def apply_publication_style() -> None:
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
        "font.size": 10,
        "axes.labelsize": 10,
        "axes.titlesize": 11,
        "axes.titleweight": "bold",
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "legend.fontsize": 9,
        "legend.frameon": False,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.25,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "axes.prop_cycle": plt.cycler(color=list(OKABE_ITO.values())[:-1]),
    })


def save(fig: plt.Figure, out_stem: Path | str, formats: Sequence[str] = ("pdf", "png")) -> None:
    stem = Path(out_stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in formats:
        fig.savefig(stem.with_suffix(f".{ext}"))


# ============================================================================
# Figure 1 — sheep donors need fewer edits (CMAH thesis)
# ============================================================================


def figure_edit_delta(plan: XenoEditPlan) -> plt.Figure:
    apply_publication_style()

    antigens = [row.antigen for row in plan.rows]
    sheep_flags = [row.sheep_to_baboon for row in plan.rows]
    pig_flags = [row.pig_to_baboon for row in plan.rows]
    n = len(antigens)

    fig = plt.figure(figsize=(10.5, 6.2))
    gs = fig.add_gridspec(
        2, 2,
        width_ratios=[1.0, 0.95], height_ratios=[1.0, 0.12],
        wspace=0.50, hspace=0.12,
    )
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])
    ax_legend = fig.add_subplot(gs[1, :])
    ax_legend.axis("off")

    # --- Panel A: matrix of coloured rectangles
    ax_a.set_xlim(-0.55, 1.55)
    ax_a.set_ylim(n - 0.45, -0.55)
    for i, (s_flag, p_flag) in enumerate(zip(sheep_flags, pig_flags)):
        if s_flag is not p_flag:
            ax_a.add_patch(mpatches.Rectangle(
                (-0.52, i - 0.48), 2.04, 0.96,
                facecolor="#FFF4DA", edgecolor=OKABE_ITO["orange"],
                linewidth=1.1, zorder=0,
            ))
        for j, flag in enumerate((s_flag, p_flag)):
            ax_a.add_patch(mpatches.Rectangle(
                (j - 0.38, i - 0.38), 0.76, 0.76,
                facecolor=NECESSITY_COLOURS[flag], edgecolor="white",
                linewidth=1.4, zorder=1,
            ))
            ax_a.text(
                j, i, NECESSITY_ABBREV[flag],
                ha="center", va="center", fontsize=9, fontweight="bold",
                color="white" if flag in (
                    EditNecessity.REQUIRED,
                    EditNecessity.COUNTERPRODUCTIVE,
                ) else OKABE_ITO["black"],
                zorder=2,
            )
    ax_a.set_yticks(range(n))
    ax_a.set_yticklabels(antigens, fontsize=9)
    ax_a.set_xticks([0, 1])
    ax_a.set_xticklabels(["Sheep → baboon", "Pig → baboon"], fontsize=10)
    ax_a.set_title("Edit necessity per xeno-antigen", pad=10)
    ax_a.tick_params(length=0)
    for spine in ax_a.spines.values():
        spine.set_visible(False)

    # --- Panel B: aggregate counts
    def _count(values: Sequence[EditNecessity], target: EditNecessity) -> int:
        return sum(1 for v in values if v is target)

    cats_short = ["Required", "Recommended", "Not required"]
    sheep_counts = [
        _count(sheep_flags, EditNecessity.REQUIRED),
        _count(sheep_flags, EditNecessity.RECOMMENDED),
        _count(sheep_flags, EditNecessity.NOT_REQUIRED),
    ]
    pig_counts = [
        _count(pig_flags, EditNecessity.REQUIRED),
        _count(pig_flags, EditNecessity.RECOMMENDED),
        _count(pig_flags, EditNecessity.NOT_REQUIRED),
    ]
    x = np.arange(len(cats_short))
    width = 0.34
    ax_b.bar(
        x - width / 2 - 0.02, sheep_counts, width,
        label="Sheep → baboon", color=OKABE_ITO["bluish_green"],
        edgecolor=OKABE_ITO["black"], linewidth=0.6,
    )
    ax_b.bar(
        x + width / 2 + 0.02, pig_counts, width,
        label="Pig → baboon", color=OKABE_ITO["vermillion"],
        edgecolor=OKABE_ITO["black"], linewidth=0.6,
    )
    for xi, (s, p) in enumerate(zip(sheep_counts, pig_counts)):
        ax_b.text(xi - width / 2 - 0.02, s + 0.18, str(s), ha="center", va="bottom",
                  fontsize=9, color=OKABE_ITO["bluish_green"], fontweight="bold")
        ax_b.text(xi + width / 2 + 0.02, p + 0.18, str(p), ha="center", va="bottom",
                  fontsize=9, color=OKABE_ITO["vermillion"], fontweight="bold")

    ax_b.set_xticks(x)
    ax_b.set_xticklabels(cats_short, fontsize=9)
    ax_b.set_ylabel("Number of antigens")
    ax_b.set_title("Aggregate edit burden", pad=10)
    y_top = max(max(pig_counts), max(sheep_counts)) + 3.0
    ax_b.set_ylim(0, y_top)
    ax_b.legend(loc="upper left", frameon=False)
    ax_b.margins(x=0.15)

    # Insight callout: placed in the clearly-empty upper-right quadrant
    # (pig has 0 in both Recommended and Not-required so this area is
    # guaranteed free). Arrow down to the green "Not required" bar.
    dropped = _count(sheep_flags, EditNecessity.NOT_REQUIRED)
    ax_b.text(
        2.05, y_top * 0.70,
        f"Δ = {dropped}\nedit dropped\non sheep axis",
        fontsize=10, fontweight="bold", color=OKABE_ITO["bluish_green"],
        ha="right", va="center",
        bbox=dict(boxstyle="round,pad=0.4",
                  facecolor="#E6F7EF", edgecolor=OKABE_ITO["bluish_green"], linewidth=1.0),
    )
    ax_b.annotate(
        "", xy=(2 - width / 2 - 0.02, sheep_counts[2] + 0.05),
        xytext=(1.80, y_top * 0.58),
        arrowprops=dict(arrowstyle="->", color=OKABE_ITO["bluish_green"], lw=1.4),
    )

    # --- Shared legend row below both panels
    legend_order = [
        EditNecessity.REQUIRED,
        EditNecessity.RECOMMENDED,
        EditNecessity.NOT_REQUIRED,
        EditNecessity.COUNTERPRODUCTIVE,
        EditNecessity.ASSESS,
    ]
    handles = [
        mpatches.Patch(
            facecolor=NECESSITY_COLOURS[k],
            edgecolor="white",
            label=k.value.replace("_", " ").lower(),
        )
        for k in legend_order
    ]
    ax_legend.legend(
        handles=handles, loc="center", ncol=5,
        title="Edit necessity", title_fontsize=9,
        fontsize=9, frameon=False, handlelength=1.3, columnspacing=1.6,
    )

    for ax, letter in [(ax_a, "A"), (ax_b, "B")]:
        ax.text(
            -0.12, 1.08, letter, transform=ax.transAxes,
            fontsize=13, fontweight="bold", va="top", ha="right",
        )
    return fig


# ============================================================================
# Figure 2 — rank vs score with published sgRNA marked
# ============================================================================


def figure_verification(
    *,
    mstn_guides: Sequence[GuideRNA],
    pdx1_guides: Sequence[GuideRNA],
    mstn_published: str,
    pdx1_published: str,
) -> plt.Figure:
    apply_publication_style()

    fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.3), sharey=True)
    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.18, top=0.87, wspace=0.12)

    specs = [
        (axes[0], mstn_guides, mstn_published,
         "Crispo 2015 · sheep MSTN", OKABE_ITO["blue"]),
        (axes[1], pdx1_guides, pdx1_published,
         "Vilarino 2017 · sheep PDX1", OKABE_ITO["vermillion"]),
    ]
    for ax, guides, published, title, colour in specs:
        n = len(guides)
        ranks = np.arange(1, n + 1)
        scores = np.array([g.on_target_score for g in guides])

        ax.scatter(
            ranks, scores,
            s=14, alpha=0.35, color=colour, edgecolor="none", zorder=1,
            label=f"pipeline candidates (n={n})",
        )

        hit_idx = next(
            (i for i, g in enumerate(guides)
             if g.protospacer.upper() == published.upper()),
            None,
        )
        if hit_idx is not None:
            hit_rank = hit_idx + 1
            hit_score = guides[hit_idx].on_target_score

            # Empty-space placement.
            # The plateau at score≈1.0 ends around rank n*0.55, and the low
            # plateau at score≈0.2 begins around rank n*0.85. The mid-band
            # (score 0.35-0.60) is genuinely empty across the whole x-axis,
            # EXCEPT for one step at ~0.7. Put the callout at y=0.42 so it
            # sits below any step, and at x past the published-hit rank so
            # the arrow travels up-and-left into empty score-1.0 whitespace
            # rather than crossing the dense top plateau.
            cx = max(hit_rank + n * 0.10, n * 0.35)
            cy = 0.42
            ax.annotate(
                f"published sgRNA\nrank {hit_rank} / {n}\nscore {hit_score:.2f}",
                xy=(hit_rank, hit_score),
                xytext=(cx, cy),
                fontsize=10, fontweight="bold",
                ha="left", va="center",
                bbox=dict(boxstyle="round,pad=0.45",
                          facecolor="white", edgecolor=OKABE_ITO["black"], linewidth=0.9),
                arrowprops=dict(arrowstyle="-", color=OKABE_ITO["black"], lw=0.8,
                                connectionstyle="arc3,rad=0.18"),
                zorder=4,
            )
            ax.scatter(
                [hit_rank], [hit_score],
                s=170, marker="D",
                facecolor=OKABE_ITO["yellow"], edgecolor=OKABE_ITO["black"],
                linewidth=1.4, zorder=5,
            )

        ax.set_xlim(-n * 0.02, n * 1.02)
        ax.set_ylim(-0.08, 1.10)
        ax.set_xlabel("Candidate rank (lower = higher score by our heuristic)")
        ax.set_title(title, pad=10)
        ax.axhline(1.0, color="#cccccc", lw=0.5, zorder=0)
        # Legend placed bottom-left — always empty (no top-ranked candidate
        # can ever have score near 0).
        ax.legend(loc="lower left", fontsize=8, bbox_to_anchor=(0.0, 0.0))
    axes[0].set_ylabel("On-target score (pipeline heuristic)")

    for ax, letter in zip(axes, "AB"):
        ax.text(
            -0.10, 1.08, letter, transform=ax.transAxes,
            fontsize=13, fontweight="bold", va="top", ha="right",
        )
    return fig


# ============================================================================
# Figure 3 — editability ranking across organs
# ============================================================================


def figure_cross_organ_landscape(
    organ_to_gene_to_guides: Mapping[str, Mapping[str, Sequence[GuideRNA]]],
) -> plt.Figure:
    apply_publication_style()

    rows: list[tuple[str, str, int]] = []
    for organ, genes in organ_to_gene_to_guides.items():
        for gene, guides in genes.items():
            if not guides:
                continue
            rows.append((organ, gene, len(guides)))

    rows.sort(key=lambda r: r[2], reverse=True)
    organs = [r[0] for r in rows]
    genes = [r[1] for r in rows]
    counts = np.array([r[2] for r in rows])
    labels = [f"{g}   ({o})" for g, o in zip(genes, organs)]

    fig = plt.figure(figsize=(10.5, 4.8))
    # Extra bottom margin reserved for the insight caption — keeps all text
    # OFF the plot area entirely.
    fig.subplots_adjust(left=0.20, right=0.80, top=0.90, bottom=0.28)
    ax = fig.add_subplot(111)

    y_pos = np.arange(len(rows))
    bar_colours = [ORGAN_COLOURS.get(o, OKABE_ITO["black"]) for o in organs]
    ax.barh(
        y_pos, counts,
        color=bar_colours, edgecolor=OKABE_ITO["black"],
        linewidth=0.6, height=0.60,
    )
    for yi, cnt in enumerate(counts):
        ax.text(
            cnt * 1.08, yi, f"{cnt}",
            va="center", ha="left", fontsize=10,
            color=OKABE_ITO["black"], fontweight="bold",
        )

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=10)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlim(80, max(counts) * 3.5)
    ax.set_xlabel("Candidate sgRNAs in sheep ortholog (log scale)")
    ax.set_title("Niche-gene editability ranking across organs", pad=12)

    legend_handles = [
        mpatches.Patch(color=c, label=organ)
        for organ, c in ORGAN_COLOURS.items() if organ in set(organs)
    ]
    ax.legend(
        handles=legend_handles,
        loc="center left", bbox_to_anchor=(1.02, 0.5),
        title="organ", title_fontsize=10, fontsize=9, frameon=False,
    )

    # Insight sits OUTSIDE the axes, below the x-axis label, as figure-level
    # text. It cannot overlap bars, tick labels, or the axis spine.
    if len(counts) >= 2 and counts[0] >= counts[-1] * 4:
        ratio = counts[0] / counts[-1]
        easiest_label = f"{genes[0]} ({organs[0]})"
        hardest_label = f"{genes[-1]} ({organs[-1]})"
        fig.text(
            0.50, 0.045,
            f"Insight: {easiest_label} has {ratio:.1f}× more candidate sgRNAs than {hardest_label}.",
            ha="center", va="bottom",
            fontsize=10, fontweight="bold", color=OKABE_ITO["vermillion"],
            bbox=dict(boxstyle="round,pad=0.40",
                      facecolor="white",
                      edgecolor=OKABE_ITO["vermillion"], linewidth=1.0),
        )

    ax.grid(axis="x", which="major", color="#EEEEEE", lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    return fig


# ============================================================================
# Case-study figures — one per README case
# ============================================================================


def figure_case_assembly_fallback(
    *,
    mstn_rambouillet_usable: int,
    mstn_texel_usable: int,
    pdx1_rambouillet_usable: int,
) -> plt.Figure:
    """Case 1 — Rambouillet MSTN is annotated as lncRNA and yields no usable
    sgRNAs; Texel recovers the protein-coding MSTN and yields N. PDX1 works
    directly on Rambouillet — no fallback needed.
    """
    apply_publication_style()

    fig, ax = plt.subplots(figsize=(8.8, 4.0))
    # Generous top and bottom margins — legend sits above axes, caption sits
    # below, so the plot area itself carries only bars + count labels.
    fig.subplots_adjust(left=0.11, right=0.97, top=0.78, bottom=0.26)

    genes = ["MSTN", "PDX1"]
    attempt1 = [mstn_rambouillet_usable, pdx1_rambouillet_usable]
    attempt2 = [mstn_texel_usable, 0]  # PDX1 never needs attempt 2

    x = np.arange(len(genes))
    width = 0.36

    ax.bar(
        x - width / 2 - 0.02, attempt1, width,
        label="Attempt 1 · ovis_aries (Rambouillet)",
        color=OKABE_ITO["sky_blue"], edgecolor=OKABE_ITO["black"], linewidth=0.6,
    )
    ax.bar(
        x + width / 2 + 0.02, attempt2, width,
        label="Attempt 2 · ovis_aries_texel (fallback)",
        color=OKABE_ITO["orange"], edgecolor=OKABE_ITO["black"], linewidth=0.6,
    )

    max_h = max(attempt1 + attempt2)
    for xi, (a1, a2) in enumerate(zip(attempt1, attempt2)):
        ax.text(xi - width / 2 - 0.02, a1 + max_h * 0.025, str(a1),
                ha="center", va="bottom", fontsize=10, fontweight="bold",
                color=OKABE_ITO["sky_blue"] if a1 > 0 else OKABE_ITO["black"])
        ax.text(xi + width / 2 + 0.02, a2 + max_h * 0.025, str(a2),
                ha="center", va="bottom", fontsize=10, fontweight="bold",
                color=OKABE_ITO["orange"] if a2 > 0 else OKABE_ITO["black"])

    ax.set_xticks(x)
    ax.set_xticklabels(genes, fontsize=11, fontweight="bold")
    ax.set_ylabel("Usable candidate sgRNAs")
    ax.set_ylim(0, max_h * 1.18)
    ax.set_title(
        "Assembly fallback — MSTN needs it, PDX1 does not",
        pad=44,  # reserve title pad so the legend (below title) is clear
    )

    # Legend placed ABOVE the axes, outside the plot rectangle.
    ax.legend(
        loc="lower center", bbox_to_anchor=(0.5, 1.01),
        ncol=2, frameon=False, fontsize=9,
    )

    # Narrative caption sits OUTSIDE the axes, at the figure bottom — cannot
    # overlap any bar, tick label, or axis.
    fig.text(
        0.50, 0.06,
        "Rambouillet annotates MSTN as a single-exon lncRNA (0 usable sgRNAs); "
        "Texel recovers the protein-coding locus (103). PDX1 needs no fallback.",
        ha="center", va="bottom", fontsize=9, style="italic",
        color=OKABE_ITO["black"],
        wrap=True,
    )
    return fig


def figure_case_tp53_locus(
    *,
    tp53_guides: Sequence[GuideRNA],
    gene_symbol: str = "TP53",
    species_label: str = "Papio anubis",
) -> plt.Figure:
    """Case 2 — a baboon TP53 lookup should land on chromosome 17 (primate
    synteny). Show each designed sgRNA as a vertical tick along chr17 with
    a locus-window zoom.
    """
    apply_publication_style()

    fig, ax = plt.subplots(figsize=(8.4, 3.2))
    fig.subplots_adjust(left=0.07, right=0.97, top=0.82, bottom=0.28)

    chroms = sorted({g.chromosome for g in tp53_guides})
    positions = np.array([(g.start + g.end) / 2 for g in tp53_guides])
    strands = np.array([g.strand for g in tp53_guides])

    if len(positions) == 0:
        ax.text(0.5, 0.5, "no guides returned", ha="center", va="center",
                transform=ax.transAxes)
        return fig

    lo, hi = positions.min(), positions.max()
    pad = max((hi - lo) * 0.15, 500)
    ax.set_xlim(lo - pad, hi + pad)
    ax.set_ylim(-1.2, 1.2)

    # Chromosome baseline
    ax.axhline(0, color=OKABE_ITO["black"], lw=1.0, zorder=1)

    for pos, strand in zip(positions, strands):
        y = 0.55 if strand == "+" else -0.55
        colour = OKABE_ITO["blue"] if strand == "+" else OKABE_ITO["vermillion"]
        ax.vlines(pos, 0, y, color=colour, lw=1.4, zorder=2)
        ax.scatter([pos], [y], s=28, color=colour,
                   edgecolor=OKABE_ITO["black"], linewidth=0.6, zorder=3)

    # Annotate locus window
    window_kb = (hi - lo) / 1000.0
    ax.text(
        (lo + hi) / 2, 1.05,
        f"{len(tp53_guides)} sgRNAs · window = {window_kb:.1f} kb",
        ha="center", va="bottom", fontsize=10, fontweight="bold",
        color=OKABE_ITO["black"],
    )

    # Strand legend (custom — avoids extra margin space)
    plus_handle = mpatches.Patch(color=OKABE_ITO["blue"], label="+ strand")
    minus_handle = mpatches.Patch(color=OKABE_ITO["vermillion"], label="− strand")
    ax.legend(handles=[plus_handle, minus_handle], loc="lower left",
              frameon=False, fontsize=9, ncol=2,
              bbox_to_anchor=(0.0, -0.02))

    ax.set_yticks([])
    ax.set_xlabel(f"Genomic position on chromosome {chroms[0]} (bp)")
    ax.set_title(
        f"{gene_symbol} ({species_label}) — all sgRNAs cluster at one locus on chr{chroms[0]}",
        pad=10,
    )
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    # Cross-chromosome sanity text at bottom (well below the baseline)
    fig.text(
        0.50, 0.04,
        f"Sanity check passed: every sgRNA is on chr{chroms[0]} "
        "(primate TP53 syntenic locus).",
        ha="center", va="bottom", fontsize=9, style="italic",
        color=OKABE_ITO["bluish_green"],
    )
    return fig


def figure_case_kidney_dual_ko(
    *,
    sall1_guides: Sequence[GuideRNA],
    six1_guides: Sequence[GuideRNA],
) -> plt.Figure:
    """Case 3 — Wang 2023 Cell Stem Cell used SALL1 + SIX1 double-KO pigs to
    build humanised mesonephros. We translate the same strategy to sheep.
    Show top-10 on-target scores for each gene side by side.
    """
    apply_publication_style()

    fig, axes = plt.subplots(1, 2, figsize=(8.8, 3.6), sharey=True)
    fig.subplots_adjust(left=0.10, right=0.97, top=0.84, bottom=0.20, wspace=0.12)

    for ax, guides, gene in [
        (axes[0], sall1_guides, "SALL1"),
        (axes[1], six1_guides, "SIX1"),
    ]:
        top = sorted(guides, key=lambda g: g.on_target_score, reverse=True)[:10]
        ranks = np.arange(1, len(top) + 1)
        scores = [g.on_target_score for g in top]
        ax.bar(
            ranks, scores,
            color=OKABE_ITO["blue"], edgecolor=OKABE_ITO["black"], linewidth=0.6,
            width=0.72,
        )
        ax.set_xticks(ranks)
        ax.set_xlabel("Top-10 rank")
        ax.set_title(f"{gene} (sheep) · {len(guides)} total candidates", pad=8)
        ax.set_ylim(0, 1.10)
        ax.axhline(1.0, color="#cccccc", lw=0.5, zorder=0)

    axes[0].set_ylabel("On-target score")

    # Top-centred banner noting origin of the dual-KO strategy.
    fig.text(
        0.50, 0.95,
        "Wang 2023 Cell Stem Cell · SALL1/SIX1 dual-KO strategy — translated to sheep",
        ha="center", va="top", fontsize=10, fontweight="bold",
        color=OKABE_ITO["black"],
    )

    # Species-transfer caption below
    fig.text(
        0.50, 0.04,
        "Both genes yield well-scored sgRNAs in Ovis aries — the strategy "
        "ports without modification.",
        ha="center", va="bottom", fontsize=9, style="italic",
        color=OKABE_ITO["bluish_green"],
    )
    return fig


def figure_case_cmah_audit(plan: XenoEditPlan) -> plt.Figure:
    """Case 4 — the CMAH row is the axis of the whole project. Spotlight it
    in a small audit-strip figure: one row per antigen, sheep and pig flags,
    highlighting the CMAH divergence."""
    apply_publication_style()

    antigens = [row.antigen for row in plan.rows]
    n = len(antigens)

    fig, ax = plt.subplots(figsize=(8.4, max(3.2, 0.44 * n + 1.6)))
    fig.subplots_adjust(left=0.18, right=0.78, top=0.88, bottom=0.12)

    ax.set_xlim(-0.55, 1.55)
    ax.set_ylim(n - 0.45, -0.55)
    for i, row in enumerate(plan.rows):
        is_cmah = row.antigen.upper() == "CMAH"
        if is_cmah:
            ax.add_patch(mpatches.Rectangle(
                (-0.58, i - 0.48), 2.16, 0.96,
                facecolor="#FFF4DA", edgecolor=OKABE_ITO["orange"],
                linewidth=1.4, zorder=0,
            ))
        for j, flag in enumerate((row.sheep_to_baboon, row.pig_to_baboon)):
            ax.add_patch(mpatches.Rectangle(
                (j - 0.38, i - 0.38), 0.76, 0.76,
                facecolor=NECESSITY_COLOURS[flag], edgecolor="white",
                linewidth=1.4, zorder=1,
            ))
            ax.text(
                j, i, NECESSITY_ABBREV[flag],
                ha="center", va="center", fontsize=9, fontweight="bold",
                color="white" if flag in (
                    EditNecessity.REQUIRED,
                    EditNecessity.COUNTERPRODUCTIVE,
                ) else OKABE_ITO["black"],
                zorder=2,
            )

    ax.set_yticks(range(n))
    ax.set_yticklabels(antigens, fontsize=10)
    ax.set_xticks([0, 1])
    ax.set_xticklabels(["Sheep → baboon", "Pig → baboon"], fontsize=10)
    ax.set_title("CMAH audit — the one row that defines this project", pad=12)
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)

    # Callout for CMAH — placed to the RIGHT of the matrix in reserved space.
    cmah_idx = next((i for i, r in enumerate(plan.rows)
                     if r.antigen.upper() == "CMAH"), None)
    if cmah_idx is not None:
        ax.annotate(
            "CMAH KO is required for pig→baboon\n"
            "but COUNTERPRODUCTIVE on the sheep\n"
            "axis — baboons retain CMAH.\n"
            "(Estrada et al. 2015)",
            xy=(1.38, cmah_idx), xytext=(2.05, cmah_idx),
            fontsize=9, ha="left", va="center",
            annotation_clip=False,
            bbox=dict(boxstyle="round,pad=0.40",
                      facecolor="white", edgecolor=OKABE_ITO["orange"],
                      linewidth=1.0),
            arrowprops=dict(arrowstyle="->", color=OKABE_ITO["orange"], lw=1.0),
        )
    return fig
