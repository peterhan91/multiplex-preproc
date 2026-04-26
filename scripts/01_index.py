#!/usr/bin/env python3
"""Inventory CODEX/CyCIF/IMC/MIBI cubes + paired H&E for the active dataset.

Builds preproc/<dataset>/index.json:
    {
      "dataset_id": ...,
      "modality": ...,
      "slides": {
        "<sample_id>": {
          "codex": {"path": ..., "shape": ..., "dtype": ..., "channels": [...], ...},
          "he":    {"path": ..., "shape": ..., "dtype": ..., ...}
        }, ...
      }
    }

The CODEX-side glob and H&E-side patterns come from the dataset config so the
same script handles Greenwald-style 1-file-per-patient OME-TIFF and
Lin2022-style WD-XXXXX-NNN per-section cubes (just by changing the YAML).
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import tifffile

from _common import (
    INDEX_JSON,
    ROOT,
    get_logger,
    load_config,
    load_index,
    patient_id_from_codex,
    save_index,
)

log = get_logger("index")


def read_ome_meta(path: Path) -> dict:
    with tifffile.TiffFile(path) as tf:
        s0 = tf.series[0]
        shape = list(s0.shape)
        axes = s0.axes
        dtype = str(s0.dtype)
        ome_xml = tf.ome_metadata or ""
        channels = re.findall(r"<Channel[^>]*Name=\"([^\"]+)\"", ome_xml)
        ps_x = ps_y = None
        m = re.search(r"PhysicalSizeX=\"([0-9.eE+-]+)\"", ome_xml)
        if m:
            ps_x = float(m.group(1))
        m = re.search(r"PhysicalSizeY=\"([0-9.eE+-]+)\"", ome_xml)
        if m:
            ps_y = float(m.group(1))
        return {
            "shape": shape,
            "axes": axes,
            "dtype": dtype,
            "channels": channels,
            "pixel_size_um": ps_x or ps_y,
        }


def read_he_meta(path: Path) -> dict:
    suffix = path.suffix.lower()
    if suffix in {".tif", ".tiff"}:
        try:
            with tifffile.TiffFile(path) as tf:
                s0 = tf.series[0]
                shape = list(s0.shape)
                axes = s0.axes
                dtype = str(s0.dtype)
                ome_xml = tf.ome_metadata or ""
                ps = None
                m = re.search(r"PhysicalSizeX=\"([0-9.eE+-]+)\"", ome_xml)
                if m:
                    ps = float(m.group(1))
                if ps is None:
                    p0 = tf.pages[0]
                    res = p0.tags.get("XResolution")
                    rx = res.value if res else None
                    unit = p0.tags.get("ResolutionUnit")
                    if rx and isinstance(rx, tuple) and rx[1]:
                        if unit and getattr(unit, "value", None) == 3:
                            ps = 10000.0 * rx[1] / rx[0]
                        else:
                            ps = 25400.0 * rx[1] / rx[0]
                return {
                    "shape": shape, "axes": axes, "dtype": dtype,
                    "pixel_size_um": ps, "format": "tiff",
                }
        except Exception as e:
            return {"format": "tiff_unreadable", "error": str(e)}
    return {"format": "unsupported", "suffix": suffix}


def discover_he(sample_id: str, he_root: Path, he_suffixes: list[str]) -> Path | None:
    if not he_root.exists():
        return None
    candidates = []
    sid_norm = sample_id.replace("_", "").lower()
    for p in he_root.rglob("*"):
        if not p.is_file():
            continue
        if p.suffix.lower() not in set(he_suffixes):
            continue
        if sid_norm in p.stem.replace("_", "").lower():
            candidates.append(p)
    if not candidates:
        return None
    candidates.sort(key=lambda p: len(p.name))
    return candidates[0]


def main() -> int:
    dc = load_config()
    raw = dc.raw_dir
    if not raw.exists():
        log.error("raw dir does not exist: %s", raw)
        return 1

    cfg_idx = dc.section("index")
    codex_glob = cfg_idx.get("codex_glob") or cfg_idx.get("cube_glob") or "*.ome.tif"
    he_subdir = cfg_idx.get("he_subdir", "he_raw")
    he_suffixes = cfg_idx.get("he_suffixes") or [".tif", ".tiff"]
    he_root = dc.preproc_dir / he_subdir

    # Start each indexing run from a clean slate so stale sids from prior runs
    # (e.g. before a config change) don't linger and cause downstream errors.
    index = {"dataset_id": dc.dataset_id, "modality": dc.modality, "slides": {}}
    slides = index["slides"]

    cube_files = sorted(raw.glob(codex_glob))
    if not cube_files:
        log.warning("no cube files matching %s in %s", codex_glob, raw)

    for codex_path in cube_files:
        sid = patient_id_from_codex(codex_path.name, dc)
        log.info("indexing CODEX %s (sample_id=%s)", codex_path.name, sid)
        codex_meta = read_ome_meta(codex_path)
        codex_meta["path"] = str(codex_path.relative_to(ROOT))
        slides.setdefault(sid, {})["codex"] = codex_meta

        he_path = discover_he(sid, he_root, he_suffixes)
        if he_path is None:
            log.info("  no matching H&E found for %s", sid)
            slides[sid]["he"] = None
        else:
            log.info("  paired H&E: %s", he_path.name)
            he_meta = read_he_meta(he_path)
            he_meta["path"] = str(he_path.relative_to(ROOT))
            slides[sid]["he"] = he_meta

    save_index(index, dc)
    log.info("wrote %s with %d slide(s)", INDEX_JSON, len(slides))
    print(json.dumps({"dataset": dc.dataset_id, "slides": list(slides.keys())}, indent=2))
    return 0


if __name__ == "__main__":
    sys.exit(main())
