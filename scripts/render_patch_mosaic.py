#!/usr/bin/env python3
"""Regenerate patch preview mosaics from a patch HDF5 file.

Default output matches the preview produced by ``04_patchify.py``:
an 8 x 8 grid of `[H&E | channel grayscale | H&E + cyan channel overlay]`
tiles using channel 0 (DAPI / nuclear reference).
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import h5py
import numpy as np

from _common import ROOT, get_logger, load_config
from mosaics import render_triptych_mosaic, resize_max_width, save_png


log = get_logger("render-mosaic")


def _parse_color(value: str) -> tuple[int, int, int]:
    parts = [p.strip() for p in value.split(",")]
    if len(parts) != 3:
        raise argparse.ArgumentTypeError("color must be R,G,B")
    try:
        rgb = tuple(int(p) for p in parts)
    except ValueError as e:
        raise argparse.ArgumentTypeError("color must contain integer R,G,B values") from e
    if any(v < 0 or v > 255 for v in rgb):
        raise argparse.ArgumentTypeError("color values must be in [0, 255]")
    return rgb


def _default_output(sample: str, channel: int, dataset_id: str) -> Path:
    return ROOT / "results" / "mosaics" / dataset_id / f"{sample}_ch{channel:02d}_preview.png"


def render_from_h5(
    h5_path: Path,
    out_path: Path,
    channel: int,
    n_patches: int,
    cols: int,
    color: tuple[int, int, int],
    alpha: float,
    lo_p: float,
    hi_p: float,
    max_width: int | None,
) -> None:
    with h5py.File(h5_path, "r") as f:
        g = f["patches"]
        he_ds = g["he"]
        codex_ds = g["codex"]
        n_total = he_ds.shape[0]
        n_channels = codex_ds.shape[1]
        if channel < 0 or channel >= n_channels:
            raise IndexError(f"channel {channel} out of range for {n_channels} channels")

        n = min(n_patches, n_total)
        he_patches = []
        channel_patches = []
        for i in range(n):
            he_patches.append(np.asarray(he_ds[i]))
            channel_patches.append(np.asarray(codex_ds[i, channel]))

        mosaic = render_triptych_mosaic(
            he_patches,
            channel_patches,
            cols=cols,
            color=color,
            alpha=alpha,
            lo_p=lo_p,
            hi_p=hi_p,
        )
        mosaic = resize_max_width(mosaic, max_width)
        save_png(mosaic, out_path)

        channel_names = list(f.attrs.get("channels", []))
        channel_name = channel_names[channel] if channel < len(channel_names) else f"ch{channel}"
        if isinstance(channel_name, bytes):
            channel_name = channel_name.decode()
        log.info(
            "wrote %s from %s: %d/%d patches, channel %d (%s), shape=%s",
            out_path,
            h5_path,
            n,
            n_total,
            channel,
            channel_name,
            mosaic.shape,
        )


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample", required=True, help="sample id, e.g. ZH1041_T1")
    ap.add_argument("--h5", type=Path, help="explicit patch HDF5 path")
    ap.add_argument("--output", type=Path, help="output PNG path")
    ap.add_argument("--channel", type=int, default=0, help="multiplex channel index")
    ap.add_argument("--n-patches", type=int, default=64, help="number of HDF5 patches to render")
    ap.add_argument("--cols", type=int, default=8, help="source patches per row")
    ap.add_argument("--color", type=_parse_color, default=(0, 255, 255), help="overlay RGB as R,G,B")
    ap.add_argument("--alpha", type=float, default=0.6, help="overlay opacity")
    ap.add_argument("--lo-p", type=float, default=1.0, help="low percentile for channel display")
    ap.add_argument("--hi-p", type=float, default=99.5, help="high percentile for channel display")
    ap.add_argument("--max-width", type=int, help="optional output width cap for shareable PNGs")
    args = ap.parse_args()

    dc = load_config()
    h5_path = args.h5 or dc.patches_dir / f"{args.sample}.h5"
    out_path = args.output or _default_output(args.sample, args.channel, dc.dataset_id)
    if not h5_path.exists():
        log.error("missing %s; run 04_patchify.py first", h5_path)
        return 1
    if args.n_patches <= 0:
        log.error("--n-patches must be positive")
        return 1
    if args.cols <= 0:
        log.error("--cols must be positive")
        return 1
    if args.max_width is not None and args.max_width <= 0:
        log.error("--max-width must be positive")
        return 1
    render_from_h5(
        h5_path=h5_path,
        out_path=out_path,
        channel=args.channel,
        n_patches=args.n_patches,
        cols=args.cols,
        color=args.color,
        alpha=args.alpha,
        lo_p=args.lo_p,
        hi_p=args.hi_p,
        max_width=args.max_width,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
