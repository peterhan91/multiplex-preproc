"""Shared helpers — now per-dataset-config-aware.

The active dataset is selected via the DATASET environment variable; default
is `greenwald-gbm-codex` (preserves the original single-dataset workflow).
Path constants are derived from the loaded dataset config so each dataset
gets its own raw/<id>/ and preproc/<id>/ subtrees.

Approach
--------
Cross-section paired multiplex IF (CODEX/CyCIF/PhenoCycler/IMC/MIBI) ↔ H&E.
The original Greenwald case shipped without a config; everything was hardcoded
under preproc/ and raw/. After the v0.5 multi-dataset refactor (2026-04-26),
that case is now described by `configs/datasets/greenwald-gbm-codex.yaml` and
its outputs live under `preproc/greenwald-gbm-codex/`. New datasets are added
by dropping a YAML next to it.
"""

from __future__ import annotations

import json
import logging
import os
import sys
from dataclasses import dataclass
from pathlib import Path

from dataset_config import (  # noqa: F401  (re-export for callers)
    DEFAULT_DATASET,
    DatasetConfig,
    active_dataset_id,
    ensure_dirs,
    load_config,
)

# Project root is the parent of scripts/. ROOT is the only thing exported as a
# module-level constant; per-dataset paths come from the loaded DatasetConfig.
ROOT = Path(__file__).resolve().parent.parent

# Back-compat module-level path constants. These are bound to the *active*
# dataset on import so existing scripts that still reference RAW/PREPROC/etc.
# work without code change. New code should prefer `load_config().<path>`.
_DC = load_config()
ensure_dirs(_DC)
RAW = _DC.raw_dir
PREPROC = _DC.preproc_dir
HE_RAW = _DC.he_raw_dir
THUMBS = _DC.thumbs_dir
REG = _DC.reg_dir
PATCHES = _DC.patches_dir
INDEX_JSON = _DC.index_json


def get_logger(name: str) -> logging.Logger:
    logging.basicConfig(
        level=os.environ.get("LOGLEVEL", "INFO"),
        format="%(asctime)s %(levelname)s %(name)s | %(message)s",
        stream=sys.stderr,
    )
    return logging.getLogger(name)


# ----------------------------------------------------------------------- device

@dataclass
class Device:
    use_cuda: bool
    torch_available: bool
    kornia_available: bool

    def __str__(self) -> str:
        return f"Device(cuda={self.use_cuda}, torch={self.torch_available}, kornia={self.kornia_available})"


def detect_device() -> Device:
    torch_ok = False
    cuda_ok = False
    kornia_ok = False
    try:
        import torch  # noqa
        torch_ok = True
        cuda_ok = bool(getattr(torch, "cuda", None) and torch.cuda.is_available())
    except Exception:
        pass
    try:
        import kornia  # noqa
        kornia_ok = True
    except Exception:
        pass
    return Device(use_cuda=cuda_ok, torch_available=torch_ok, kornia_available=kornia_ok)


# ----------------------------------------------------------------------- index

def load_index(dc: DatasetConfig | None = None) -> dict:
    fp = (dc or _DC).index_json
    if fp.exists():
        return json.loads(fp.read_text())
    return {"slides": {}}


def save_index(index: dict, dc: DatasetConfig | None = None) -> None:
    (dc or _DC).index_json.write_text(json.dumps(index, indent=2, sort_keys=True))


# --------------------------------------------------------------- patient match

def patient_id_from_codex(codex_filename: str, dc: DatasetConfig | None = None) -> str:
    """Map a CODEX/IMC/MIBI filename → patient/sample id.

    Strip rules come from the dataset config:
      index.patient_id_strip_suffixes: list[str]    (e.g. ['_B','_A'])
      index.patient_id_strip_version_re: regex      (e.g. '^v\\d+$')

    Greenwald default: ZH1041_T1_B.ome.tif → ZH1041_T1; MGH258_v6.ome.tif → MGH258.
    """
    import re
    dc = dc or _DC
    name = Path(codex_filename).name
    # Strip extension. ".ome.tif" is a double suffix so Path.stem alone leaves
    # ".ome"; handle that case explicitly. Falls through to Path.stem for other
    # formats (.mcd, .tiff, .qptiff, etc.).
    if ".ome.tif" in name:
        stem = name.split(".ome.tif")[0]
    else:
        stem = Path(name).stem
    parts = stem.split("_")
    cfg_idx = dc.section("index")
    strip_suffixes = cfg_idx.get("patient_id_strip_suffixes") or []
    strip_ver_re = cfg_idx.get("patient_id_strip_version_re")
    if parts and parts[-1] in set(strip_suffixes):
        parts = parts[:-1]
    if parts and strip_ver_re and re.fullmatch(strip_ver_re, parts[-1]):
        parts = parts[:-1]
    return "_".join(parts)


def patient_id_from_he(he_filename: str, dc: DatasetConfig | None = None) -> str:
    """Map an H&E filename → sample id space (mirror of patient_id_from_codex).

    Strips the common H&E decoration suffixes after the standard CODEX-side
    strips have been applied.
    """
    name = Path(he_filename).stem
    for tag in ("_H&E", "_HE", "_he", "_h_e", "-HE", "-he"):
        if name.endswith(tag):
            name = name[: -len(tag)]
            break
    return name
