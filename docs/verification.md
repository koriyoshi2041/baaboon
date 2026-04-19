# Verification strategy

This document specifies *how we know the pipeline is correct*. The goal is
explicit, auditable, tiered checks rather than "it ran". Every Tier is
implemented as code under `tests/` and runs in CI.

## Tier 0 — Static catalog invariants

Offline checks that every YAML file conforms to structural expectations,
every curated gene carries a citation, and the project's central
scientific claim (CMAH delta for sheep→baboon) is encoded correctly.

Tests: `tests/test_catalog.py`, `tests/test_xeno_delta.py`.
Runtime: < 1 s. No network.

## Tier 1 — Golden-standard recapitulation

For each published edit in `data/golden_standards.yaml`, we run the pipeline
on the same locus and species, collect the top-N scored sgRNAs, and assert
that the published protospacer is present with rank ≤ N (default N=10).

Published sources currently listed:

| Study | Species | Gene | Strategy |
| --- | --- | --- | --- |
| Vilarino et al. 2017 *Sci Rep* | *Ovis aries* | PDX1 | dual sgRNA, exon 1 |
| Crispo et al. 2015 *PLOS One* | *Ovis aries* | MSTN | single sgRNA |
| Niu et al. 2017 *Science* | *Sus scrofa* | PERV-pol | 1 sgRNA × 59 loci |

The sequences themselves are transcribed from each paper's supplementary
materials. Until that transcription is complete (`status: pending_curation`
in the YAML), the test `test_vilarino_pdx1_guides_recover_published` is
`xfail` with a clear message so unattended curation is visible.

Tests: `tests/test_vilarino_pdx1.py`.
Runtime: Ensembl-bound, < 30 s per gene with cache warm.

## Tier 2 — Cross-tool consistency

Off by default; opt-in via `pytest -m integration`. Runs the pipeline
alongside CRISPOR and CHOPCHOP on the same target and asserts:

- top-20 Jaccard index ≥ 0.4
- Spearman rank correlation of scores ≥ 0.6

CRISPOR and CHOPCHOP are external tools not bundled here; the harness
invokes their CLIs if `CRISPOR_PATH` / `CHOPCHOP_PATH` env vars point to
working installations.

## Tier 3 — Biological sanity

Lightweight assertions that the plan output respects known biology:

1. **Phylogenetic identity ordering.** For any gene where all three species
   annotate an ortholog, human↔baboon protein identity must exceed
   sheep↔baboon identity. If this inverts, the ortholog resolver is
   broken.
2. **CMAH demotion.** In the Stage 2 table, `sheep_to_baboon_edit` for
   CMAH must be `NOT_REQUIRED`, reflecting Estrada 2015.
3. **TP53 competence uplift.** If/when a chimerism simulator is added,
   TP53-KO plans should score higher than wild-type plans.
4. **Cell adhesion awareness.** Stage 1b must always emit at least one
   plan from the `cell_adhesion` class, reflecting the Lin 2024 addition
   to the barrier landscape.

Tests: `tests/test_xeno_delta.py` for items 2 and 4. Items 1 and 3 are
planned follow-ups.

## Golden-curation checklist

To move a `golden_standards.yaml` entry from `pending_curation` to
`curated`:

1. Download the paper's supplementary material (PDF or Excel).
2. Locate the exact sgRNA table. Note the column labels; distinguish
   **target / protospacer** from **full sgRNA + scaffold**.
3. Extract the 20 nt protospacer in 5'→3' orientation. Confirm the PAM
   is NGG adjacent (not included in the field).
4. Paste into the YAML as `protospacer:`.
5. Change `status: pending_curation` → `status: curated`.
6. Run `pytest -m golden -m network`. Expect the matching test to pass.
7. Commit with a message citing the paper DOI.

## How to report a verification failure

If a Tier 1 or Tier 3 check fails unexpectedly:

1. Record the pipeline version (`papovis --version`).
2. Record the Ensembl release (pinned in `ASSEMBLY_BY_SPECIES`).
3. Include the raw pipeline output for the failing target.
4. File an issue tagged `verification`.

Failures of Tier 3 are especially informative — they suggest either a
curation mistake or a genuine scientific update to the literature.
