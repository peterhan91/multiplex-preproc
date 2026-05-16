"""Patch mosaic rendering helpers.

These functions are intentionally small and data-format agnostic: callers pass
already extracted H&E RGB patches plus one multiplex channel.  The same helpers
are used by the pipeline preview and the standalone regeneration CLI so QC
images do not drift across code paths.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
from PIL import Image


RGBColor = tuple[int, int, int]


def to_u8(arr: np.ndarray, lo_p: float = 1.0, hi_p: float = 99.5) -> np.ndarray:
    """Percentile-normalize one image plane to uint8."""
    a = arr.astype(np.float32)
    lo, hi = np.percentile(a, [lo_p, hi_p])
    if hi <= lo:
        hi = lo + 1.0
    return (np.clip((a - lo) / (hi - lo), 0, 1) * 255).astype(np.uint8)


def overlay_channel(
    he_rgb: np.ndarray,
    ch_u8: np.ndarray,
    color: RGBColor = (0, 255, 255),
    alpha: float = 0.6,
) -> np.ndarray:
    """Additively overlay a uint8 protein/nuclear channel on an H&E patch."""
    color_arr = np.zeros_like(he_rgb)
    color_arr[..., 0] = (ch_u8 * (color[0] / 255.0)).astype(np.uint8)
    color_arr[..., 1] = (ch_u8 * (color[1] / 255.0)).astype(np.uint8)
    color_arr[..., 2] = (ch_u8 * (color[2] / 255.0)).astype(np.uint8)
    blend = np.clip(
        he_rgb.astype(np.uint16) + (color_arr.astype(np.uint16) * alpha).astype(np.uint16),
        0,
        255,
    )
    return blend.astype(np.uint8)


def triptych_tile(
    he_rgb: np.ndarray,
    channel_plane: np.ndarray,
    color: RGBColor = (0, 255, 255),
    alpha: float = 0.6,
    lo_p: float = 1.0,
    hi_p: float = 99.5,
) -> np.ndarray:
    """Render one `[H&E | channel grayscale | H&E + channel overlay]` tile."""
    ch_u8 = to_u8(channel_plane, lo_p=lo_p, hi_p=hi_p)
    ch_rgb = np.stack([ch_u8] * 3, axis=-1)
    overlay = overlay_channel(he_rgb, ch_u8, color=color, alpha=alpha)
    return np.hstack([he_rgb, ch_rgb, overlay])


def render_triptych_mosaic(
    he_patches: Iterable[np.ndarray],
    channel_patches: Iterable[np.ndarray],
    cols: int = 8,
    color: RGBColor = (0, 255, 255),
    alpha: float = 0.6,
    lo_p: float = 1.0,
    hi_p: float = 99.5,
) -> np.ndarray:
    """Render a grid of triptych tiles.

    `cols` is the number of source patches per row.  The rendered image has
    three visual columns per source patch because each patch is a triptych.
    """
    tiles = [
        triptych_tile(he, ch, color=color, alpha=alpha, lo_p=lo_p, hi_p=hi_p)
        for he, ch in zip(he_patches, channel_patches)
    ]
    if not tiles:
        raise ValueError("cannot render an empty patch mosaic")

    rows = []
    for r0 in range(0, len(tiles), cols):
        rows.append(np.hstack(tiles[r0 : r0 + cols]))
    return np.vstack(rows)


def save_png(arr: np.ndarray, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    Image.fromarray(arr).save(out_path)


def resize_max_width(arr: np.ndarray, max_width: int | None) -> np.ndarray:
    """Downsample an RGB array to `max_width`, preserving aspect ratio."""
    if max_width is None or max_width <= 0 or arr.shape[1] <= max_width:
        return arr
    height = max(1, int(arr.shape[0] * (max_width / arr.shape[1])))
    im = Image.fromarray(arr)
    return np.asarray(im.resize((max_width, height), Image.Resampling.LANCZOS))
