"""Ensembl REST client with retry, on-disk caching, and typed responses.

The client targets `rest.ensembl.org` release 115 (September 2025), which is
the first release to include a freshly annotated Panubis1.0 (olive baboon).
Callers should not rely on implicit Ensembl version; the client pins
`Accept: application/json` and includes the release number in the cache key.
"""

from __future__ import annotations

import hashlib
import json
import logging
import os
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import requests
from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential,
)

from papovis.models import ASSEMBLY_BY_SPECIES, SequenceRegion, Species

logger = logging.getLogger(__name__)

DEFAULT_SERVER = "https://rest.ensembl.org"
DEFAULT_RELEASE = 115
DEFAULT_CACHE_DIR = Path(
    os.environ.get("PAPOVIS_CACHE_DIR", str(Path.home() / ".cache" / "papovis"))
)

# Ensembl REST returns human-species_name strings; map back from our Species literal.
_SPECIES_TO_ENSEMBL_NAME: dict[Species, str] = {
    "Ovis aries": "ovis_aries",
    "Papio anubis": "papio_anubis",
    "Homo sapiens": "homo_sapiens",
    "Sus scrofa": "sus_scrofa",
}

# Fallback assembly aliases to use when the primary species name at Ensembl
# release 115 has incomplete annotation for the queried gene. Concrete
# observed example: sheep MSTN is annotated as lncRNA on the Rambouillet
# reference but as protein_coding on the Texel reference.
_SPECIES_FALLBACKS: dict[Species, list[str]] = {
    "Ovis aries": ["ovis_aries_texel"],
}


class EnsemblError(RuntimeError):
    """Non-retryable Ensembl REST error (4xx other than rate-limit)."""


class EnsemblClient:
    """Typed, cached Ensembl REST wrapper."""

    def __init__(
        self,
        *,
        server: str = DEFAULT_SERVER,
        release: int = DEFAULT_RELEASE,
        cache_dir: Path = DEFAULT_CACHE_DIR,
        session: requests.Session | None = None,
        timeout: float = 30.0,
    ) -> None:
        self._server = server.rstrip("/")
        self._release = release
        self._cache_dir = Path(cache_dir)
        self._cache_dir.mkdir(parents=True, exist_ok=True)
        self._session = session or requests.Session()
        self._timeout = timeout

    # ------------------------------------------------------------------ HTTP

    def _cache_key(self, path: str, params: Mapping[str, Any] | None) -> Path:
        payload = json.dumps(
            {"release": self._release, "path": path, "params": dict(params or {})},
            sort_keys=True,
        )
        digest = hashlib.sha256(payload.encode()).hexdigest()
        return self._cache_dir / f"{digest}.json"

    @retry(
        retry=retry_if_exception_type((requests.ConnectionError, requests.Timeout)),
        stop=stop_after_attempt(4),
        wait=wait_exponential(multiplier=1.0, min=1.0, max=16.0),
        reraise=True,
    )
    def _http_get(self, path: str, params: Mapping[str, Any] | None) -> Any:
        url = f"{self._server}{path}"
        response = self._session.get(
            url,
            params=dict(params or {}),
            headers={"Accept": "application/json"},
            timeout=self._timeout,
        )
        if response.status_code == 429:
            # Signal a retry via a Timeout (which tenacity will retry on).
            raise requests.Timeout(f"Ensembl rate limit on {url}")
        if 400 <= response.status_code < 500:
            raise EnsemblError(
                f"Ensembl REST {response.status_code} for {url}: {response.text[:200]}"
            )
        response.raise_for_status()
        return response.json()

    def _get(self, path: str, params: Mapping[str, Any] | None = None) -> Any:
        cache_path = self._cache_key(path, params)
        if cache_path.exists():
            try:
                return json.loads(cache_path.read_text())
            except json.JSONDecodeError:
                logger.warning("Corrupt cache at %s, re-fetching", cache_path)
                cache_path.unlink(missing_ok=True)
        payload = self._http_get(path, params)
        cache_path.write_text(json.dumps(payload))
        return payload

    # -------------------------------------------------------------- Endpoints

    def lookup_gene_by_symbol(self, species: Species, gene_symbol: str) -> dict[str, Any]:
        """Resolve a gene symbol to an Ensembl gene record.

        If the primary assembly returns no protein_coding transcripts (a
        known Ensembl annotation gap, e.g. sheep MSTN on Rambouillet), the
        client transparently retries against the registered fallback
        assembly names. The returned record adds a
        ``_papovis_resolved_species`` key naming the assembly that was
        actually used, for downstream logging.
        """
        candidates = [_SPECIES_TO_ENSEMBL_NAME[species]] + _SPECIES_FALLBACKS.get(species, [])
        last_error: EnsemblError | None = None
        for ensembl_name in candidates:
            try:
                payload = self._get(f"/lookup/symbol/{ensembl_name}/{gene_symbol}", {"expand": "1"})
            except EnsemblError as exc:
                last_error = exc
                continue
            if not isinstance(payload, dict) or "id" not in payload:
                continue
            if _has_protein_coding_transcript(payload):
                payload["_papovis_resolved_species"] = ensembl_name
                if ensembl_name != candidates[0]:
                    logger.info(
                        "Gene %s resolved via fallback assembly %s for %s.",
                        gene_symbol,
                        ensembl_name,
                        species,
                    )
                return payload
            # Remember this as a weak hit in case no fallback yields a coding transcript.
            last_error = EnsemblError(
                f"{gene_symbol} on {ensembl_name} has no protein_coding transcript."
            )
        if last_error is not None:
            raise last_error
        raise EnsemblError(
            f"Ensembl did not return a gene record for {gene_symbol} in {species}."
        )

    def fetch_sequence(
        self,
        species: Species,
        chromosome: str,
        start: int,
        end: int,
        strand: int = 1,
        ensembl_name_override: str | None = None,
    ) -> SequenceRegion:
        """Fetch genomic DNA for a 1-based inclusive region.

        `ensembl_name_override` lets the caller pin a specific assembly
        alias (e.g. ``ovis_aries_texel``) so that coordinates resolved from
        one assembly are not accidentally queried against another.
        """
        ensembl_name = ensembl_name_override or _SPECIES_TO_ENSEMBL_NAME[species]
        region = f"{chromosome}:{start}..{end}:{strand}"
        path = f"/sequence/region/{ensembl_name}/{region}"
        payload = self._get(path)
        sequence = payload["seq"] if isinstance(payload, dict) else payload[0]["seq"]
        return SequenceRegion(
            species=species,
            chromosome=chromosome,
            start=start,
            end=end,
            strand="+" if strand == 1 else "-",
            sequence=sequence.upper(),
        )

    def fetch_exon_sequences(
        self,
        species: Species,
        gene_symbol: str,
        max_exons: int = 3,
    ) -> list[SequenceRegion]:
        """Fetch forward-strand sequence of the first N exons of the canonical
        protein-coding transcript.

        Coordinates are resolved through `lookup_gene_by_symbol`, which may
        route to a fallback assembly alias when the primary assembly has
        incomplete annotation. The fetched sequence uses the SAME assembly
        that resolved the gene, so coordinates and sequences are always
        consistent.
        """
        gene = self.lookup_gene_by_symbol(species, gene_symbol)
        resolved = gene.get(
            "_papovis_resolved_species", _SPECIES_TO_ENSEMBL_NAME[species]
        )
        canonical = _pick_canonical_transcript(gene)
        exons = canonical.get("Exon", [])[:max_exons]
        regions: list[SequenceRegion] = []
        for exon in exons:
            regions.append(
                self.fetch_sequence(
                    species=species,
                    chromosome=str(exon["seq_region_name"]),
                    start=int(exon["start"]),
                    end=int(exon["end"]),
                    strand=1,
                    ensembl_name_override=resolved,
                )
            )
        return regions

    def fetch_orthologs(
        self,
        source_species: Species,
        gene_symbol: str,
        target_species: Species,
    ) -> list[dict[str, Any]]:
        """Query Ensembl Compara for orthologs of a gene between two species."""
        source_name = _SPECIES_TO_ENSEMBL_NAME[source_species]
        target_name = _SPECIES_TO_ENSEMBL_NAME[target_species]
        path = f"/homology/symbol/{source_name}/{gene_symbol}"
        payload = self._get(
            path,
            {
                "target_species": target_name,
                "type": "orthologues",
                "format": "condensed",
            },
        )
        data = payload.get("data", []) if isinstance(payload, dict) else []
        if not data:
            return []
        return data[0].get("homologies", [])


def _has_protein_coding_transcript(gene_record: dict[str, Any]) -> bool:
    """True if the gene record contains at least one protein_coding transcript."""
    for transcript in gene_record.get("Transcript", []):
        if transcript.get("biotype") == "protein_coding":
            return True
    return False


def _pick_canonical_transcript(gene_record: dict[str, Any]) -> dict[str, Any]:
    """Pick a canonical transcript from an expanded gene lookup payload.

    Preference order (protein_coding strongly preferred because non-coding
    canonical-flagged transcripts are common Ensembl artefacts — e.g. sheep
    MSTN at release 115 flags a 1-exon lncRNA as canonical while the
    multi-exon protein-coding transcript is unflagged):

      1. Transcript with biotype == "protein_coding" and is_canonical == 1
      2. Any transcript with is_canonical == 1
      3. Longest protein_coding transcript by CDS length
      4. Longest transcript overall
    """
    transcripts = gene_record.get("Transcript", [])
    if not transcripts:
        raise EnsemblError("Gene record has no transcripts.")

    def is_protein_coding(t: dict[str, Any]) -> bool:
        return t.get("biotype") == "protein_coding"

    canonical_pc = [
        t for t in transcripts if t.get("is_canonical") == 1 and is_protein_coding(t)
    ]
    if canonical_pc:
        return canonical_pc[0]

    canonicals = [t for t in transcripts if t.get("is_canonical") == 1]
    protein_coding = [t for t in transcripts if is_protein_coding(t)]
    if protein_coding:
        # Prefer longest CDS when multiple protein-coding transcripts exist.
        def cds_length(t: dict[str, Any]) -> int:
            translation = t.get("Translation") or {}
            return int(translation.get("length", 0))

        return max(protein_coding, key=cds_length)

    if canonicals:
        return canonicals[0]

    def total_length(t: dict[str, Any]) -> int:
        return int(t.get("end", 0)) - int(t.get("start", 0))

    return max(transcripts, key=total_length)
