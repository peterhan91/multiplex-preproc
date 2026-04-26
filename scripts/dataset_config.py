"""Per-dataset YAML config loader.

A dataset config lives at `configs/datasets/<dataset_id>.yaml` and is the
single source of truth for: which raw archives to extract, which CODEX/IMC/MIBI
files to index, which channel is nuclear-reference, which registration
strategy to use, where the panel comes from, and patchify parameters.

The active dataset is selected via the DATASET environment variable; pipelines
default to the constant ``DEFAULT_DATASET`` below if it isn't set, so existing
single-dataset workflows keep working.
"""
from __future__ import annotations

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml

DEFAULT_DATASET = "greenwald-gbm-codex"

ROOT = Path(__file__).resolve().parent.parent
CONFIGS_DIR = ROOT / "configs" / "datasets"


@dataclass
class DatasetConfig:
    raw: dict
    cfg: dict

    # convenience accessors so call sites read like attributes, not dict lookups
    @property
    def dataset_id(self) -> str:
        return self.cfg["dataset_id"]

    @property
    def modality(self) -> str:
        return self.cfg["modality"]

    @property
    def raw_dir(self) -> Path:
        sub = self.cfg.get("paths", {}).get("raw_subdir") or self.dataset_id
        return ROOT / "raw" / sub

    @property
    def preproc_dir(self) -> Path:
        sub = self.cfg.get("paths", {}).get("preproc_subdir") or self.dataset_id
        return ROOT / "preproc" / sub

    @property
    def index_json(self) -> Path:
        return self.preproc_dir / "index.json"

    @property
    def thumbs_dir(self) -> Path:
        return self.preproc_dir / "thumbs"

    @property
    def reg_dir(self) -> Path:
        return self.preproc_dir / "registration"

    @property
    def patches_dir(self) -> Path:
        return self.preproc_dir / "patches"

    @property
    def panels_dir(self) -> Path:
        return self.preproc_dir / "panels"

    @property
    def logs_dir(self) -> Path:
        return self.preproc_dir / "logs"

    @property
    def he_raw_dir(self) -> Path:
        sub = self.cfg.get("index", {}).get("he_subdir", "he_raw")
        return self.preproc_dir / sub

    def section(self, key: str) -> dict:
        return self.cfg.get(key) or {}


def active_dataset_id() -> str:
    return os.environ.get("DATASET", DEFAULT_DATASET)


def load_config(dataset_id: str | None = None) -> DatasetConfig:
    did = dataset_id or active_dataset_id()
    fp = CONFIGS_DIR / f"{did}.yaml"
    if not fp.exists():
        raise SystemExit(
            f"dataset config not found: {fp}\n"
            f"Available: {sorted(p.stem for p in CONFIGS_DIR.glob('*.yaml') if not p.stem.startswith('_'))}"
        )
    raw = yaml.safe_load(fp.read_text())
    cfg = _validate(raw, fp)
    return DatasetConfig(raw=raw, cfg=cfg)


def _validate(raw: dict, fp: Path) -> dict:
    required_top = {"dataset_id", "modality"}
    missing = required_top - set(raw or {})
    if missing:
        raise SystemExit(f"{fp}: missing required keys: {sorted(missing)}")
    valid_modalities = {"codex", "phenocycler", "cycif", "orion", "imc", "mibi"}
    if raw["modality"] not in valid_modalities:
        raise SystemExit(
            f"{fp}: invalid modality {raw['modality']!r} — must be one of {sorted(valid_modalities)}"
        )
    return raw


def ensure_dirs(dc: DatasetConfig) -> None:
    """Create the per-dataset preproc subdirs the pipeline writes into."""
    for d in (dc.preproc_dir, dc.thumbs_dir, dc.reg_dir, dc.patches_dir,
              dc.panels_dir, dc.logs_dir, dc.he_raw_dir):
        d.mkdir(parents=True, exist_ok=True)
