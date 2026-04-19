"""Load the curated YAML catalogs of niche genes, barrier genes, xeno antigens,
and golden-standard sgRNAs. Everything the pipeline reads from disk flows
through this module.
"""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path
from typing import Any

import yaml

DATA_DIR = Path(__file__).resolve().parent.parent / "data"


def _load_yaml(name: str) -> dict[str, Any]:
    path = DATA_DIR / name
    if not path.exists():
        raise FileNotFoundError(f"Missing data file: {path}")
    with path.open("r", encoding="utf-8") as fh:
        payload = yaml.safe_load(fh)
    if not isinstance(payload, dict):
        raise ValueError(f"Expected a mapping in {path}, got {type(payload).__name__}")
    return payload


@lru_cache(maxsize=1)
def niche_catalog() -> dict[str, Any]:
    return _load_yaml("niche_genes.yaml")


@lru_cache(maxsize=1)
def barrier_catalog() -> dict[str, Any]:
    return _load_yaml("barrier_genes.yaml")


@lru_cache(maxsize=1)
def xeno_catalog() -> dict[str, Any]:
    return _load_yaml("xeno_antigens.yaml")


@lru_cache(maxsize=1)
def golden_catalog() -> dict[str, Any]:
    return _load_yaml("golden_standards.yaml")


def available_organs() -> list[str]:
    return sorted(niche_catalog().get("organs", {}).keys())


def niche_genes_for_organ(organ: str) -> list[dict[str, Any]]:
    organs = niche_catalog().get("organs", {})
    if organ not in organs:
        raise KeyError(
            f"Unknown organ {organ!r}. Known organs: {', '.join(available_organs())}"
        )
    return list(organs[organ].get("niche_genes", []))


def barrier_genes(class_name: str | None = None) -> list[dict[str, Any]]:
    classes = barrier_catalog().get("classes", {})
    if class_name is None:
        rows: list[dict[str, Any]] = []
        for cls_name, cls in classes.items():
            for gene in cls.get("genes", []):
                rows.append({**gene, "barrier_class": cls_name})
        return rows
    if class_name not in classes:
        raise KeyError(f"Unknown barrier class {class_name!r}")
    return [
        {**gene, "barrier_class": class_name}
        for gene in classes[class_name].get("genes", [])
    ]


def xeno_antigens() -> dict[str, list[dict[str, Any]]]:
    catalog = xeno_catalog()
    return {
        "glycan_antigens": list(catalog.get("glycan_antigens", [])),
        "complement_regulators": list(catalog.get("complement_regulators", [])),
        "coagulation_thrombosis": list(catalog.get("coagulation_thrombosis", [])),
    }


def golden_entries() -> list[dict[str, Any]]:
    return list(golden_catalog().get("entries", []))
