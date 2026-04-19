"""Shared pytest fixtures.

Network-requiring tests must be marked with `@pytest.mark.network`.
By default those are deselected (see pyproject.toml) unless CI opts in with
`-m network`.
"""

from __future__ import annotations

import pytest


@pytest.fixture
def synthetic_pdx1_like_sequence() -> str:
    """A fabricated DNA region containing multiple SpCas9 targets for unit tests.

    The sequence is deterministic, deliberately contains several NGG PAMs, a
    homopolymer run, and a polyT-terminator substring, so that each scorer
    branch exercises at least one test case.
    """
    # 150 bp — enough for a handful of 20+3 protospacer+PAM hits on both strands.
    return (
        "ATGCCGACCTTGCAGGTCAAGACTGGCGAACAGCTCAAGGCAGGCATCGAT"  # 50
        "CGGATTCCAAAAAAGTGCTACCGAGAGCGGTATTCGGCCAATGGCGATGCC"  # 100
        "TGGCAGTTTAGCCCAGCTGACTTAGTGGCCTGGATCCCGAAGGTACCGTAC"  # 150
    )


@pytest.fixture
def short_sequence_no_pam() -> str:
    return "AAAAAAAAAAAAAAAAAAAA"


@pytest.fixture
def min_pam_sequence() -> str:
    """Exactly one NGG PAM at the end, with a clean 20-mer protospacer."""
    return "ACGTACGTACGTACGTACGT" + "AGG"
