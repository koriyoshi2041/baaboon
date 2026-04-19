"""Publication-quality figure generation for papovis results.

Style follows top-tier biology journal conventions:
  - Single column width = 89 mm (≈3.5 in); double column = 183 mm (≈7.2 in)
  - Sans-serif typography (DejaVu Sans → Helvetica → Arial fallback order)
  - Strict Okabe–Ito colour-blind-safe palette (Okabe & Ito 2008)
  - Despined axes, no gridlines, no 3-D effects
  - Bold panel letters (A, B, C) in the top-left
  - Vector PDF + 300 dpi PNG outputs
  - Insight-first: every panel has one thing to say

Each function takes prepared data objects (not raw pipeline state) so the
figures can be regenerated from cached intermediates and exercised by unit
tests without hitting the network.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

from papovis.models import EditNecessity, GuideRNA, XenoEditPlan

# ---------------------------------------------------------------- Okabe-Ito
# The eight canonical Okabe & Ito colours. Distinguishable under all common
# forms of colour-blindness (deuteranopia, protanopia, tritanopia).
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

# Semantic colour map for EditNecessity — strict subset of Okabe-Ito.
NECESSITY_COLOURS: dict[EditNecessity, str] = {
    EditNecessity.REQUIRED: OKABE_ITO["vermillion"],
    EditNecessity.RECOMMENDED: OKABE_ITO["orange"],
    EditNecessity.ASSESS: OKABE_ITO["sky_blue"],
    EditNecessity.NOT_REQUIRED: OKABE_ITO["bluish_green"],
    EditNecessity.COUNTERPRODUCTIVE: OKABE_ITO["reddish_purple"],
}

# Short tokens for in-cell labels (≤3 chars, readable at small size).
NECESSITY_ABBREV: dict[EditNecessity, str] = {
    EditNecessity.REQUIRED: "REQ",
    EditNecessity.RECOMMENDED: "rec",
    EditNecessity.ASSESS: "?",
    EditNecessity.NOT_REQUIRED: "—",
    EditNecessity.COUNTERPRODUCTIVE: "✗",
}

# Organ → Okabe-Ito colour (used for Figure 3 bar colouring).
ORGAN_COLOURS: dict[str, str] = {
    "pancreas": OKABE_ITO["vermillion"],
    "kidney": OKABE_ITO["blue"],
    "thymus": OKABE_ITO["bluish_green"],
    "liver": OKABE_ITO["orange"],
    "heart": OKABE_ITO["reddish_purple"],
}


def apply_publication_style() -> None:
    """Install the project-wide matplotlib rcParams."""
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
        "font.size": 9,
        "axes.labelsize": 9,
        "axes.titlesize": 10,
        "axes.titleweight": "bold",
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "legend.fontsize": 8,
        "legend.frameon": False,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.fonttype": "none",
        "axes.prop_cycle": plt.cycler(color=list(OKABE_ITO.values())[:-1]),
    })


def save(fig: plt.Figure, out_stem: Path | str, formats: Sequence[str] = ("pdf", "png")) -> None:
    """Save a figure to multiple formats using the publication style."""
    stem = Path(out_stem)
    stem.parent.mkdir(parents=True, exist_ok=True)
    for ext in formats:
        fig.savefig(stem.with_suffix(f".{ext}"))


# ============================================================================
# Figure 1 — "sheep donors need fewer edits" (CMAH thesis)
# ============================================================================


def figure_edit_delta(plan: XenoEditPlan) -> plt.Figure:
    """Two panels.

    **A.** Per-antigen side-by-side heatmap. Rows = antigens, two columns =
    (sheep axis, pig axis). Cell colour = edit-necessity category. The CMAH
    row is visually emphasised because its necessity flips between axes —
    the project's signature finding.

    **B.** Aggregate count bar chart. Sheep vs pig, for each necessity
    category. The Δ is annotated.
    """
    apply_publication_style()

    antigens = [row.antigen for row in plan.rows]
    sheep_flags = [row.sheep_to_baboon for row in plan.rows]
    pig_flags = [row.pig_to_baboon for row in plan.rows]
    n = len(antigens)

    fig = plt.figure(figsize=(7.2, max(3.8, 0.28 * n + 1.6)))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.15, 1.0], wspace=0.55)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])

    # --- Panel A: matrix of coloured rectangles
    ax_a.set_xlim(-0.5, 1.5)
    ax_a.set_ylim(n - 0.5, -0.5)
    for i, (s_flag, p_flag) in enumerate(zip(sheep_flags, pig_flags)):
        # Highlight the row where sheep vs pig diverge — the CMAH insight.
        if s_flag is not p_flag:
            ax_a.add_patch(mpatches.Rectangle(
                (-0.48, i - 0.48), 1.96, 0.96,
                facecolor="#FFF3D6", edgecolor="#E69F00", linewidth=1.2, zorder=0,
            ))
        for j, flag in enumerate((s_flag, p_flag)):
            colour = NECESSITY_COLOURS[flag]
            ax_a.add_patch(mpatches.Rectangle(
                (j - 0.40, i - 0.40), 0.80, 0.80,
                facecolor=colour, edgecolor="white", linewidth=1.5, zorder=1,
            ))
            ax_a.text(
                j, i, NECESSITY_ABBREV[flag],
                ha="center", va="center", fontsize=8, fontweight="bold",
                color="white" if flag in (
                    EditNecessity.REQUIRED,
                    EditNecessity.COUNTERPRODUCTIVE,
                    EditNecessity.RECOMMENDED,
                ) else OKABE_ITO["black"],
                zorder=2,
            )
    ax_a.set_yticks(range(n))
    ax_a.set_yticklabels(antigens, fontsize=8)
    ax_a.set_xticks([0, 1])
    ax_a.set_xticklabels(["Sheep → baboon", "Pig → baboon"])
    ax_a.set_title("Edit necessity per xeno-antigen")
    ax_a.tick_params(length=0)
    for spine in ax_a.spines.values():
        spine.set_visible(False)

    # Legend above / beside.
    legend_order = [
        EditNecessity.REQUIRED,
        EditNecessity.RECOMMENDED,
        EditNecessity.ASSESS,
        EditNecessity.NOT_REQUIRED,
        EditNecessity.COUNTERPRODUCTIVE,
    ]
    handles = [
        mpatches.Patch(
            facecolor=NECESSITY_COLOURS[k],
            edgecolor="white",
            label=k.value.replace("_", " ").lower(),
        )
        for k in legend_order
    ]
    ax_a.legend(
        handles=handles,
        loc="upper left",
        bbox_to_anchor=(1.02, 1.0),
        borderaxespad=0.2,
        fontsize=7,
        title="necessity",
        title_fontsize=7,
    )

    # --- Panel B: aggregate bars
    def _count(values: Sequence[EditNecessity], target: EditNecessity) -> int:
        return sum(1 for v in values if v is target)

    cats = ["REQUIRED", "RECOMMENDED", "NOT\nREQUIRED"]
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
    x = np.arange(len(cats))
    width = 0.36
    ax_b.bar(
        x - width / 2, sheep_counts, width,
        label="Sheep → baboon", color=OKABE_ITO["bluish_green"],
        edgecolor=OKABE_ITO["black"], linewidth=0.6,
    )
    ax_b.bar(
        x + width / 2, pig_counts, width,
        label="Pig → baboon", color=OKABE_ITO["vermillion"],
        edgecolor=OKABE_ITO["black"], linewidth=0.6,
    )
    # Value labels on top of bars
    for xi, (s, p) in enumerate(zip(sheep_counts, pig_counts)):
        ax_b.text(xi - width / 2, s + 0.12, str(s), ha="center", va="bottom",
                  fontsize=7, color=OKABE_ITO["bluish_green"], fontweight="bold")
        ax_b.text(xi + width / 2, p + 0.12, str(p), ha="center", va="bottom",
                  fontsize=7, color=OKABE_ITO["vermillion"], fontweight="bold")

    ax_b.set_xticks(x)
    ax_b.set_xticklabels(cats)
    ax_b.set_ylabel("Number of antigens")
    ax_b.set_title("Aggregate edit burden")
    ax_b.set_ylim(0, max(max(pig_counts), max(sheep_counts)) + 1.6)
    ax_b.legend(loc="upper right")

    dropped = _count(sheep_flags, EditNecessity.NOT_REQUIRED)
    # Big Δ annotation
    ax_b.annotate(
        f"Δ = {dropped}\nedit dropped\non sheep axis",
        xy=(2 - width / 2, sheep_counts[2] + 0.05),
        xytext=(1.05, max(pig_counts) * 0.55),
        fontsize=8, fontweight="bold", color=OKABE_ITO["bluish_green"],
        arrowprops=dict(arrowstyle="->", color=OKABE_ITO["bluish_green"], lw=1.2),
    )

    # Panel labels
    for ax, letter in [(ax_a, "A"), (ax_b, "B")]:
        ax.text(
            -0.18, 1.06, letter, transform=ax.transAxes,
            fontsize=12, fontweight="bold", va="top", ha="right",
        )
    return fig


# ============================================================================
# Figure 2 — rank vs score view of pipeline output, marking published sgRNAs
# ============================================================================


def figure_verification(
    *,
    mstn_guides: Sequence[GuideRNA],
    pdx1_guides: Sequence[GuideRNA],
    mstn_published: str,
    pdx1_published: str,
) -> plt.Figure:
    """Two panels, one per published experiment.

    For each, every candidate sgRNA is plotted as a dot (rank on x, score
    on y). The published protospacer is marked with a large diamond and a
    callout box stating rank and score. The insight: how far a transparent
    heuristic scorer and the original published tool disagree on ranking.
    """
    apply_publication_style()

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.2), sharey=True)

    specs = [
        (axes[0], mstn_guides, mstn_published,
         "Crispo 2015 · MSTN · sheep", OKABE_ITO["blue"]),
        (axes[1], pdx1_guides, pdx1_published,
         "Vilarino 2017 · PDX1 · sheep", OKABE_ITO["vermillion"]),
    ]
    for ax, guides, published, title, colour in specs:
        ranks = np.arange(1, len(guides) + 1)
        scores = np.array([g.on_target_score for g in guides])

        # All candidates as translucent dots
        ax.scatter(
            ranks, scores,
            s=10, alpha=0.35, color=colour,
            edgecolor="none", zorder=1,
            label=f"pipeline candidates (n={len(guides)})",
        )

        # Published guide as large diamond
        hit_idx = next(
            (i for i, g in enumerate(guides)
             if g.protospacer.upper() == published.upper()),
            None,
        )
        if hit_idx is not None:
            hit_rank = hit_idx + 1
            hit_score = guides[hit_idx].on_target_score
            ax.scatter(
                [hit_rank], [hit_score],
                s=140, marker="D",
                facecolor=OKABE_ITO["yellow"],
                edgecolor=OKABE_ITO["black"],
                linewidth=1.2, zorder=3,
                label="published sgRNA",
            )
            # Callout
            ax.annotate(
                f"rank {hit_rank} / {len(guides)}\nscore {hit_score:.2f}",
                xy=(hit_rank, hit_score),
                xytext=(min(hit_rank + len(guides) * 0.12, len(guides) * 0.6),
                        max(hit_score - 0.35, 0.05)),
                fontsize=8, fontweight="bold",
                bbox=dict(boxstyle="round,pad=0.3",
                          facecolor="white", edgecolor=OKABE_ITO["black"], linewidth=0.8),
                arrowprops=dict(arrowstyle="->", color=OKABE_ITO["black"], lw=0.8),
                zorder=4,
            )

        ax.set_xlim(0, len(guides) + 1)
        ax.set_ylim(-0.05, 1.08)
        ax.set_xlabel("Candidate rank (our scorer)")
        ax.set_title(title)
        ax.legend(loc="lower left", fontsize=7)
    axes[0].set_ylabel("On-target score")

    # Panel letters
    for ax, letter in zip(axes, "AB"):
        ax.text(
            -0.12, 1.07, letter, transform=ax.transAxes,
            fontsize=12, fontweight="bold", va="top", ha="right",
        )
    return fig


# ============================================================================
# Figure 3 — cross-organ editability ranking (REDESIGN: horizontal bars)
# ============================================================================


def figure_cross_organ_landscape(
    organ_to_gene_to_guides: Mapping[str, Mapping[str, Sequence[GuideRNA]]],
) -> plt.Figure:
    """Horizontal bar chart of candidate-sgRNA count per (organ, niche-gene)
    pair, sorted descending. Colour encodes organ. The log x-axis handles
    the ~10× range naturally and the SIX1-vs-HHEX contrast is annotated.
    """
    apply_publication_style()

    rows: list[tuple[str, str, int, float]] = []  # (organ, gene, count, mean_score)
    for organ, genes in organ_to_gene_to_guides.items():
        for gene, guides in genes.items():
            if not guides:
                continue
            top10 = sorted((g.on_target_score for g in guides), reverse=True)[:10]
            rows.append((organ, gene, len(guides), float(np.mean(top10))))

    rows.sort(key=lambda r: r[2], reverse=True)
    organs = [r[0] for r in rows]
    genes = [r[1] for r in rows]
    counts = np.array([r[2] for r in rows])
    scores = np.array([r[3] for r in rows])
    labels = [f"{g}  ·  {o}" for g, o in zip(genes, organs)]

    fig, ax = plt.subplots(figsize=(7.0, 3.4))

    y_pos = np.arange(len(rows))
    bar_colours = [ORGAN_COLOURS.get(o, OKABE_ITO["black"]) for o in organs]
    ax.barh(
        y_pos, counts,
        color=bar_colours,
        edgecolor=OKABE_ITO["black"],
        linewidth=0.6, height=0.62,
    )
    # In-bar count labels
    for yi, (cnt, score) in enumerate(zip(counts, scores)):
        ax.text(
            cnt * 1.05, yi, f"{cnt}  (score {score:.2f})",
            va="center", ha="left", fontsize=8, color=OKABE_ITO["black"],
        )

    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xscale("log")
    ax.set_xlim(80, max(counts) * 2.6)
    ax.set_xlabel("Candidate sgRNAs in sheep ortholog (log scale)")
    ax.set_title("Niche-gene editability ranking across organs")

    # Legend for organ colours
    legend_handles = [
        mpatches.Patch(color=c, label=organ)
        for organ, c in ORGAN_COLOURS.items() if organ in set(organs)
    ]
    ax.legend(
        handles=legend_handles,
        loc="lower right", title="organ",
        title_fontsize=8, fontsize=7,
    )

    # Insight annotation
    if counts[0] >= counts[-1] * 4 and len(counts) >= 2:
        ratio = counts[0] / counts[-1]
        ax.annotate(
            f"{ratio:.1f}× more candidate sgRNAs\nthan the scarcest niche target",
            xy=(counts[0], 0), xytext=(counts[0] * 1.1, len(rows) * 0.25),
            fontsize=8, fontweight="bold", color=OKABE_ITO["vermillion"],
            arrowprops=dict(arrowstyle="->", color=OKABE_ITO["vermillion"], lw=1.2),
        )

    ax.grid(axis="x", which="major", color="#E5E7EB", lw=0.6, zorder=0)
    ax.set_axisbelow(True)
    fig.tight_layout()
    return fig
