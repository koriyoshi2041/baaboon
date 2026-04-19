"""Publication-quality figure generation for papovis results.

Style follows top-tier biology journal conventions:
  - Single column width = 89 mm (≈3.5 in); double column = 183 mm (≈7.2 in)
  - Sans-serif typography (Helvetica toArial toDejaVu Sans fallback)
  - Okabe–Ito colourblind-safe palette
  - Despined axes, no gridlines, no 3-D effects
  - Bold panel letters (A, B, C) in the top-left
  - Vector PDF + 300 dpi PNG outputs

Each function takes prepared data objects (not raw pipeline state) so the
figures can be regenerated from cached intermediates and exercised by unit
tests without hitting the network.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from pathlib import Path
from string import ascii_uppercase

import matplotlib as mpl
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np

from papovis.models import EditNecessity, GuideRNA, XenoEditPlan

# Okabe-Ito colourblind-safe palette (Okabe & Ito 2008).
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

# Semantic colour map for EditNecessity.
NECESSITY_COLOURS: dict[EditNecessity, str] = {
    EditNecessity.REQUIRED: OKABE_ITO["vermillion"],
    EditNecessity.RECOMMENDED: OKABE_ITO["orange"],
    EditNecessity.ASSESS: OKABE_ITO["sky_blue"],
    EditNecessity.NOT_REQUIRED: OKABE_ITO["bluish_green"],
    EditNecessity.COUNTERPRODUCTIVE: "#6A3D9A",  # purple tail; not in Okabe-Ito
}


def apply_publication_style() -> None:
    """Install the project-wide matplotlib rcParams."""
    mpl.rcParams.update({
        # DejaVu Sans first so Unicode arrows render; Helvetica/Arial still
        # preferred when all required glyphs are present.
        "font.family": "sans-serif",
        "font.sans-serif": ["DejaVu Sans", "Helvetica", "Arial"],
        "font.size": 8,
        "axes.labelsize": 8,
        "axes.titlesize": 9,
        "axes.titleweight": "bold",
        "axes.linewidth": 0.8,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "xtick.major.width": 0.8,
        "ytick.major.width": 0.8,
        "legend.fontsize": 7,
        "legend.frameon": False,
        "figure.dpi": 150,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "pdf.fonttype": 42,  # keep text as text in PDF
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


# --------------------------------------------------------------------------
# Figure 1 — the CMAH advantage (two-panel, project thesis)
# --------------------------------------------------------------------------


def figure_edit_delta(plan: XenoEditPlan) -> plt.Figure:
    """Two-panel figure summarising the Stage-2 sheep-vs-pig edit burden.

    Panel A: heatmap of edit necessity per antigen × donor species axis.
    Panel B: grouped bar chart aggregating required/recommended/dropped edits.
    """
    apply_publication_style()

    antigens = [row.antigen for row in plan.rows]
    necessity_matrix = np.zeros((len(antigens), 2), dtype=int)
    necessity_order = [
        EditNecessity.COUNTERPRODUCTIVE,
        EditNecessity.NOT_REQUIRED,
        EditNecessity.ASSESS,
        EditNecessity.RECOMMENDED,
        EditNecessity.REQUIRED,
    ]
    necessity_index = {n: i for i, n in enumerate(necessity_order)}
    for i, row in enumerate(plan.rows):
        necessity_matrix[i, 0] = necessity_index[row.sheep_to_baboon]
        necessity_matrix[i, 1] = necessity_index[row.pig_to_baboon]

    # Build colour mapping
    cmap = mpl.colors.ListedColormap(
        [NECESSITY_COLOURS[n] for n in necessity_order]
    )

    fig = plt.figure(figsize=(7.2, 4.2))
    gs = fig.add_gridspec(1, 2, width_ratios=[1.2, 1.0], wspace=0.4)
    ax_a = fig.add_subplot(gs[0, 0])
    ax_b = fig.add_subplot(gs[0, 1])

    # --- Panel A: heatmap
    im = ax_a.imshow(
        necessity_matrix,
        aspect="auto",
        cmap=cmap,
        vmin=-0.5,
        vmax=len(necessity_order) - 0.5,
    )
    ax_a.set_yticks(range(len(antigens)))
    ax_a.set_yticklabels(antigens)
    ax_a.set_xticks([0, 1])
    ax_a.set_xticklabels(["Sheep  → baboon", "Pig  → baboon"])
    ax_a.set_xlabel("Donor axis")
    ax_a.set_title("Edit necessity per xeno-antigen")
    ax_a.tick_params(length=0)
    for spine in ax_a.spines.values():
        spine.set_visible(False)
    for i in range(len(antigens)):
        for j in range(2):
            ax_a.text(
                j, i,
                necessity_order[necessity_matrix[i, j]].value[:3],
                ha="center", va="center",
                color="white" if necessity_matrix[i, j] >= 3 else "black",
                fontsize=6,
                fontweight="bold",
            )

    # Legend for panel A
    legend_handles = [
        mpatches.Patch(color=NECESSITY_COLOURS[n], label=n.value.replace("_", " ").title())
        for n in necessity_order
    ]
    ax_a.legend(
        handles=legend_handles,
        loc="upper left",
        bbox_to_anchor=(1.0, 1.0),
        borderaxespad=0.2,
        fontsize=6,
    )

    # --- Panel B: aggregate counts
    def _count(values: Sequence[EditNecessity], target: EditNecessity) -> int:
        return sum(1 for v in values if v is target)

    sheep_values = [row.sheep_to_baboon for row in plan.rows]
    pig_values = [row.pig_to_baboon for row in plan.rows]
    categories = ["REQUIRED", "RECOMMENDED", "NOT\nREQUIRED"]
    sheep_counts = [
        _count(sheep_values, EditNecessity.REQUIRED),
        _count(sheep_values, EditNecessity.RECOMMENDED),
        _count(sheep_values, EditNecessity.NOT_REQUIRED),
    ]
    pig_counts = [
        _count(pig_values, EditNecessity.REQUIRED),
        _count(pig_values, EditNecessity.RECOMMENDED),
        _count(pig_values, EditNecessity.NOT_REQUIRED),
    ]
    x = np.arange(len(categories))
    width = 0.36
    ax_b.bar(
        x - width / 2, sheep_counts, width,
        label="Sheep  → baboon", color=OKABE_ITO["bluish_green"], edgecolor="black", linewidth=0.5,
    )
    ax_b.bar(
        x + width / 2, pig_counts, width,
        label="Pig  → baboon", color=OKABE_ITO["vermillion"], edgecolor="black", linewidth=0.5,
    )
    ax_b.set_xticks(x)
    ax_b.set_xticklabels(categories)
    ax_b.set_ylabel("Number of antigens")
    ax_b.set_title("Aggregate edit burden")
    ax_b.legend(loc="upper right")
    dropped = _count(sheep_values, EditNecessity.NOT_REQUIRED)
    ax_b.annotate(
        f"Δ = {dropped} edit(s)\ndropped on\nsheep axis",
        xy=(2, sheep_counts[2]),
        xytext=(1.3, max(pig_counts) * 0.6),
        fontsize=7,
        arrowprops=dict(arrowstyle="->", color="black", lw=0.6),
    )

    # Panel labels
    for ax, letter in [(ax_a, "A"), (ax_b, "B")]:
        ax.text(
            -0.15, 1.05, letter, transform=ax.transAxes,
            fontsize=11, fontweight="bold", va="top", ha="right",
        )
    return fig


# --------------------------------------------------------------------------
# Figure 2 — Live Tier-1 verification
# --------------------------------------------------------------------------


def figure_verification(
    *,
    mstn_guides: Sequence[GuideRNA],
    pdx1_guides: Sequence[GuideRNA],
    mstn_published: str,
    pdx1_published: str,
) -> plt.Figure:
    """Score histograms for MSTN and PDX1 runs, with the published sgRNAs
    highlighted as vertical annotations."""
    apply_publication_style()

    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0), sharey=True)

    for ax, guides, published, title, colour in [
        (axes[0], mstn_guides, mstn_published,
         "Crispo 2015 MSTN (sheep)", OKABE_ITO["blue"]),
        (axes[1], pdx1_guides, pdx1_published,
         "Vilarino 2017 PDX1 (sheep)", OKABE_ITO["vermillion"]),
    ]:
        scores = np.array([g.on_target_score for g in guides])
        ax.hist(scores, bins=np.linspace(0, 1, 21), color=colour,
                edgecolor="black", linewidth=0.4, alpha=0.85)
        hit = next((i + 1 for i, g in enumerate(guides)
                    if g.protospacer.upper() == published.upper()), None)
        if hit is not None:
            hit_score = guides[hit - 1].on_target_score
            ax.axvline(hit_score, color="black", linestyle="--", lw=0.8)
            ax.text(
                hit_score, ax.get_ylim()[1] * 0.95,
                f"  published\n  sgRNA\n  (rank {hit}/{len(guides)})",
                fontsize=6, va="top", ha="left",
            )
        ax.set_xlim(0, 1.05)
        ax.set_xlabel("On-target score (pipeline)")
        ax.set_title(title)
    axes[0].set_ylabel("Candidate sgRNAs")

    for ax, letter in zip(axes, "AB"):
        ax.text(
            -0.12, 1.06, letter, transform=ax.transAxes,
            fontsize=11, fontweight="bold", va="top", ha="right",
        )
    return fig


# --------------------------------------------------------------------------
# Figure 3 — Cross-organ niche editability landscape
# --------------------------------------------------------------------------


def figure_cross_organ_landscape(
    organ_to_gene_to_guides: Mapping[str, Mapping[str, Sequence[GuideRNA]]],
) -> plt.Figure:
    """Dual-encoded heatmap: rows = organs, columns = niche genes. Cell colour
    encodes top-10 average on-target score; cell text shows number of
    candidate guides produced by the pipeline. NaN cells indicate
    (organ, gene) pairs not curated by the catalog."""
    apply_publication_style()
    organs = list(organ_to_gene_to_guides.keys())
    all_genes = sorted({
        gene for inner in organ_to_gene_to_guides.values() for gene in inner
    })

    count_matrix = np.full((len(organs), len(all_genes)), np.nan)
    score_matrix = np.full((len(organs), len(all_genes)), np.nan)
    for i, organ in enumerate(organs):
        inner = organ_to_gene_to_guides[organ]
        for j, gene in enumerate(all_genes):
            guides = inner.get(gene)
            if not guides:
                continue
            count_matrix[i, j] = len(guides)
            top10 = sorted((g.on_target_score for g in guides), reverse=True)[:10]
            score_matrix[i, j] = float(np.mean(top10))

    fig, ax = plt.subplots(figsize=(5.5, 3.5))
    cmap = mpl.cm.viridis.copy()
    cmap.set_bad(color="#EEEEEE")
    im = ax.imshow(
        np.ma.masked_invalid(score_matrix),
        aspect="auto",
        cmap=cmap,
        vmin=0.4,
        vmax=1.0,
    )
    for i in range(len(organs)):
        for j in range(len(all_genes)):
            n = count_matrix[i, j]
            if np.isnan(n):
                continue
            ax.text(
                j, i, f"{int(n)}",
                ha="center", va="center",
                color="white" if score_matrix[i, j] >= 0.7 else "black",
                fontsize=7,
                fontweight="bold",
            )
    ax.set_xticks(range(len(all_genes)))
    ax.set_xticklabels(all_genes, rotation=45, ha="right")
    ax.set_yticks(range(len(organs)))
    ax.set_yticklabels(organs)
    ax.set_title("Niche-vacancy sgRNA landscape across organs")
    ax.tick_params(length=0)
    for spine in ax.spines.values():
        spine.set_visible(False)
    cbar = fig.colorbar(im, ax=ax, shrink=0.75, pad=0.02)
    cbar.set_label("Mean top-10 on-target score", fontsize=7)
    cbar.ax.tick_params(labelsize=6)
    fig.tight_layout()
    return fig
