#!/usr/bin/env python3
"""Per-channel QC panels: for each selected CODEX channel, render
[H&E | protein heatmap | H&E + cyan-protein overlay] across a few tissue
patches. Saves one PNG per channel to preproc/patches/<sid>_panels/.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np
from PIL import Image

from _common import get_logger, load_config
from panels import GREENWALD_PANEL

log = get_logger("panels")


def _to_u8(arr: np.ndarray, lo_p: float = 1.0, hi_p: float = 99.5) -> np.ndarray:
    a = arr.astype(np.float32)
    lo, hi = np.percentile(a, [lo_p, hi_p])
    if hi <= lo:
        hi = lo + 1.0
    return (np.clip((a - lo) / (hi - lo), 0, 1) * 255).astype(np.uint8)


def _overlay(he_rgb: np.ndarray, ch_u8: np.ndarray, color=(0, 255, 255), alpha: float = 0.6) -> np.ndarray:
    color_arr = np.zeros_like(he_rgb)
    color_arr[..., 0] = (ch_u8 * (color[0] / 255.0)).astype(np.uint8)
    color_arr[..., 1] = (ch_u8 * (color[1] / 255.0)).astype(np.uint8)
    color_arr[..., 2] = (ch_u8 * (color[2] / 255.0)).astype(np.uint8)
    blend = np.clip(he_rgb.astype(np.uint16) + (color_arr.astype(np.uint16) * alpha).astype(np.uint16), 0, 255)
    return blend.astype(np.uint8)


def panel_for_channel(
    he: np.ndarray,         # (N, P, P, 3) uint8
    codex: np.ndarray,      # (N, C, P, P) uint8
    chan_idx: int,
    chan_name: str,
    patch_indices: list[int],
    cols: int = 4,
) -> np.ndarray:
    rows = []
    for r0 in range(0, len(patch_indices), cols):
        row = []
        for i in patch_indices[r0 : r0 + cols]:
            he_p = he[i]
            ch = codex[i, chan_idx]
            ch_u8 = _to_u8(ch)
            ch_rgb = np.stack([ch_u8] * 3, axis=-1)
            ovl = _overlay(he_p, ch_u8)
            tile = np.hstack([he_p, ch_rgb, ovl])
            row.append(tile)
        if row:
            rows.append(np.hstack(row))
    return np.vstack(rows)


def _label_band(width: int, text: str, height: int = 40) -> np.ndarray:
    from PIL import Image as _Im, ImageDraw as _Dr
    im = _Im.new("RGB", (width, height), (24, 24, 24))
    dr = _Dr.Draw(im)
    dr.text((8, 10), text, fill=(240, 240, 240))
    return np.asarray(im)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True)
    ap.add_argument(
        "--channels",
        default="all",
        help="comma-separated channel indices, or 'all'",
    )
    ap.add_argument("--n_patches", type=int, default=8, help="patches per channel panel")
    ap.add_argument("--cols", type=int, default=4)
    ap.add_argument(
        "--panel",
        default="greenwald43",
        choices=["greenwald43", "fluorophore"],
        help="marker-name source: paper-verified Greenwald 43-ch panel (Table S3), or raw OME fluorophore labels",
    )
    args = ap.parse_args()

    dc = load_config()
    h5_path = dc.patches_dir / f"{args.sample}.h5"
    if not h5_path.exists():
        log.error("missing %s — run 04_patchify.py first", h5_path)
        return 1
    out_dir = dc.patches_dir / f"{args.sample}_panels"
    out_dir.mkdir(parents=True, exist_ok=True)

    with h5py.File(h5_path, "r") as f:
        g = f["patches"]
        N = g["he"].shape[0]
        C = g["codex"].shape[1]
        if args.panel == "greenwald43" and C == len(GREENWALD_PANEL):
            chan_names = list(GREENWALD_PANEL)
            log.info("using paper-verified Greenwald panel (Table S3, %d channels)", C)
        else:
            chan_names = list(f.attrs.get("channels", [f"ch{i}" for i in range(C)]))
            log.info("using OME-XML channel names (%s)", args.panel)
        if args.channels == "all":
            chan_list = list(range(C))
        else:
            chan_list = [int(c) for c in args.channels.split(",")]
        # pick patches with highest mean nuclear (DAPI) signal so panels show real tissue
        dapi_means = g["codex"][:, 0].reshape(N, -1).mean(axis=1)
        top = np.argsort(-dapi_means)[: max(args.n_patches, args.cols)]
        patch_indices = sorted(top.tolist())[: args.n_patches]
        log.info("sample %s: N=%d C=%d, panel patches: %s", args.sample, N, C, patch_indices)

        for ci in chan_list:
            if ci < 0 or ci >= C:
                log.warning("skip out-of-range channel %d", ci); continue
            name = chan_names[ci] if ci < len(chan_names) else f"ch{ci}"
            log.info("  channel %d (%s)", ci, name)
            he = g["he"][patch_indices]
            codex_sub = g["codex"][patch_indices][:, ci : ci + 1]
            # build a (N_sub, C=1, P, P) view; reshape to use the same panel function
            full_codex = np.zeros((len(patch_indices), C, codex_sub.shape[2], codex_sub.shape[3]), dtype=np.uint8)
            full_codex[:, ci] = codex_sub[:, 0]
            panel = panel_for_channel(
                he, full_codex, ci, name,
                patch_indices=list(range(len(patch_indices))),
                cols=args.cols,
            )
            label = _label_band(panel.shape[1], f"ch{ci:02d}  {name}   [H&E | protein | overlay]")
            out_img = np.vstack([label, panel])
            fp = out_dir / f"ch{ci:02d}_{name.replace(' ', '_')}.png"
            Image.fromarray(out_img).save(fp)
            log.info("    → %s", fp)

    log.info("done. panels in %s", out_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
