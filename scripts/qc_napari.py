#!/usr/bin/env python3
"""Interactive QC: load H&E + a CODEX channel as napari layers with the
registration affine baked in, so panning/zooming reveals true cell-level
alignment.

Usage
-----
    python scripts/qc_napari.py --sample ZH1041_T1 --channel DAPI
    python scripts/qc_napari.py --sample ZH1041_T1 --channel 0   # by index
    python scripts/qc_napari.py --sample ZH1041_T1 --level 1     # pyramid lvl

Requires: napari (`pip install 'napari[all]'`).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np
import tifffile

from _common import ROOT, get_logger, load_config, load_index

log = get_logger("qc_napari")


def _resolve_channel(channels: list[str], spec: str) -> int:
    if spec.isdigit():
        return int(spec)
    s = spec.lower()
    for i, c in enumerate(channels or []):
        if c and s in c.lower():
            return i
    raise SystemExit(f"channel '{spec}' not found in {channels}")


def _pyramid_level(tf: tifffile.TiffFile, level: int) -> np.ndarray:
    s0 = tf.series[0]
    if not s0.levels or level >= len(s0.levels):
        return s0.asarray()
    return s0.levels[level].asarray()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument("--channel", default="DAPI", help="CODEX channel name or index")
    ap.add_argument("--level", type=int, default=1, help="pyramid level (0=full res)")
    args = ap.parse_args()

    try:
        import napari
    except ImportError:
        log.error("napari not installed. pip install 'napari[all]'")
        return 1

    dc = load_config()
    index = load_index(dc)
    slide = index["slides"].get(args.sample)
    if not slide:
        log.error("sample %s not in index", args.sample)
        return 1
    reg_path = dc.reg_dir / f"{args.sample}.json"
    if not reg_path.exists():
        log.error("no registration for %s — run 03_register.py first", args.sample)
        return 1
    reg = json.loads(reg_path.read_text())
    if not reg.get("ok"):
        log.error("registration not ok for %s", args.sample)
        return 1

    codex_path = (ROOT / slide["codex"]["path"]).resolve()
    he_path = (ROOT / slide["he"]["path"]).resolve()
    ch_idx = _resolve_channel(slide["codex"].get("channels") or [], args.channel)
    log.info("loading CODEX ch%d (%s) and H&E at level %d", ch_idx,
             (slide["codex"].get("channels") or [None])[ch_idx], args.level)

    with tifffile.TiffFile(codex_path) as tf:
        codex_arr = _pyramid_level(tf, args.level)
        s0 = tf.series[0]
        codex_full_h = s0.shape[-2]
        codex_lvl_h = codex_arr.shape[-2]
    with tifffile.TiffFile(he_path) as tf:
        he_arr = _pyramid_level(tf, args.level)
        s0 = tf.series[0]
        he_full_h = s0.shape[-2]
        he_lvl_h = he_arr.shape[-2]

    # CODEX axes handling
    if codex_arr.ndim == 3:
        c_axis = int(np.argmin(codex_arr.shape))
        codex_ch = np.take(codex_arr, ch_idx, axis=c_axis)
    else:
        codex_ch = codex_arr

    # affine_thumb_codex_to_he is in thumbnail coords. Convert to chosen pyramid
    # level by composing scale matrices: T_lvl = S(he_lvl/he_thumb) @ T_thumb @ S(codex_thumb/codex_lvl).
    A_thumb = np.array(reg["affine_thumb_codex_to_he"], dtype=np.float64)
    s_he = he_lvl_h / reg["thumb_shape_he"][0]
    s_co_inv = reg["thumb_shape_codex"][0] / codex_lvl_h
    S_he = np.diag([s_he, s_he, 1.0])
    S_co = np.diag([s_co_inv, s_co_inv, 1.0])
    T_lvl = S_he @ A_thumb @ S_co  # CODEX_lvl → H&E_lvl

    # napari uses (y, x); affine for napari layers is on world coords.
    # The CODEX layer is in CODEX pixel space; we set its affine so it lands
    # in H&E pixel space (the H&E layer has identity affine).
    affine_yx = np.eye(3)
    affine_yx[0, 0] = T_lvl[1, 1]
    affine_yx[0, 1] = T_lvl[1, 0]
    affine_yx[0, 2] = T_lvl[1, 2]
    affine_yx[1, 0] = T_lvl[0, 1]
    affine_yx[1, 1] = T_lvl[0, 0]
    affine_yx[1, 2] = T_lvl[0, 2]

    viewer = napari.Viewer(title=f"{args.sample} — H&E + CODEX[{ch_idx}]")
    viewer.add_image(he_arr, name="H&E", rgb=he_arr.ndim == 3 and he_arr.shape[-1] in (3, 4))
    viewer.add_image(
        codex_ch,
        name=f"CODEX ch{ch_idx}",
        affine=affine_yx,
        colormap="green",
        blending="additive",
        contrast_limits=[float(np.percentile(codex_ch, 1)), float(np.percentile(codex_ch, 99.5))],
    )
    napari.run()
    return 0


if __name__ == "__main__":
    sys.exit(main())
