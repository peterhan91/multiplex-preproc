#!/usr/bin/env python3
"""Write low-res grayscale thumbnails for each registered pair (config-driven).

For the cube side: pick a pyramid level near `registration.thumb_target_downsample`
and slice the *reference channel* selected by `reference_channel.match_any`
(falling back to `fallback_index`). Different modalities use different ref
channels — DAPI for CODEX/CyCIF, Ir191/Ir193 for IMC, dsDNA for MIBI.

For the H&E: same pyramid level, RGB → 1-V (HSV) so hematoxylin-rich nuclei
become bright (matching DAPI brightness convention).

Outputs:
  preproc/<dataset>/thumbs/<sid>__codex.png  (uint8 grayscale)
  preproc/<dataset>/thumbs/<sid>__he.png     (uint8 grayscale)
  preproc/<dataset>/thumbs/<sid>__he_rgb.png (uint8 RGB, QC-only)
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import tifffile
from PIL import Image

from _common import ROOT, get_logger, load_config, load_index
from readers import resolve_reference_channel_index

log = get_logger("thumbs")


def _pick_level(tf: tifffile.TiffFile, target_ds: int) -> int:
    s0 = tf.series[0]
    base_y = s0.shape[-2]
    target_y = max(1, base_y // target_ds)
    best, best_diff = 0, abs(s0.levels[0].shape[-2] - target_y)
    for i, lvl in enumerate(s0.levels):
        d = abs(lvl.shape[-2] - target_y)
        if d < best_diff:
            best, best_diff = i, d
    return best


def _u8(arr: np.ndarray) -> np.ndarray:
    arr = arr.astype(np.float32)
    lo, hi = np.percentile(arr, [1, 99])
    if hi <= lo:
        hi = lo + 1.0
    arr = np.clip((arr - lo) / (hi - lo), 0, 1)
    return (arr * 255).astype(np.uint8)


def thumb_codex(codex_path: Path, ref_idx: int, target_ds: int) -> np.ndarray:
    """OME-TIFF path — pyramid level near target downsample, then ref channel."""
    with tifffile.TiffFile(codex_path) as tf:
        s0 = tf.series[0]
        if not s0.levels:
            arr = s0.asarray()
        else:
            lvl = _pick_level(tf, target_ds)
            log.info("  cube level %d (shape=%s)", lvl, s0.levels[lvl].shape)
            arr = s0.levels[lvl].asarray()

    if arr.ndim == 3:
        c_axis = int(np.argmin(arr.shape))
        arr = np.take(arr, ref_idx, axis=c_axis)
    elif arr.ndim == 4:
        z_axis = int(np.argmin(arr.shape))
        arr = arr.max(axis=z_axis)
        c_axis = int(np.argmin(arr.shape))
        arr = np.take(arr, ref_idx, axis=c_axis)
    return _u8(arr)


def thumb_via_reader(codex_path: Path, ref_idx: int, target_ds: int, modality: str) -> np.ndarray:
    """Modality-aware fallback for non-OME-TIFF cubes (IMC .mcd, MIBI dir).

    Opens via the reader factory, takes the ref channel, downsamples by
    `target_ds`. The reader loads the full cube into RAM (small for IMC/MIBI
    ROIs which are typically <1k × 1k px).
    """
    from readers import open_reader
    r = open_reader(codex_path, modality, preload=True)
    try:
        full = r._mem[ref_idx]  # (H, W) float32
    finally:
        r.close()
    if target_ds > 1:
        from PIL import Image as _Im
        nh = max(1, full.shape[0] // target_ds)
        nw = max(1, full.shape[1] // target_ds)
        u8 = _u8(full)
        thumb = np.asarray(_Im.fromarray(u8).resize((nw, nh), _Im.BILINEAR))
        log.info("  reader-based thumb (full %s → %s)", full.shape, thumb.shape)
        return thumb
    return _u8(full)


def thumb_he(he_path: Path, target_ds: int, method: str = "invert_value") -> tuple[np.ndarray, np.ndarray | None]:
    with tifffile.TiffFile(he_path) as tf:
        s0 = tf.series[0]
        if not s0.levels:
            arr = s0.asarray()
        else:
            lvl = _pick_level(tf, target_ds)
            log.info("  H&E level %d (shape=%s)", lvl, s0.levels[lvl].shape)
            arr = s0.levels[lvl].asarray()

    if arr.ndim == 2:
        return _u8(arr), None
    rgb = arr if arr.shape[-1] == 3 else np.moveaxis(arr, 0, -1)[..., :3]
    rgb_f = rgb.astype(np.float32) / 255.0
    if method == "invert_value":
        v = rgb_f.max(axis=-1)
        gray = 1.0 - v
    elif method == "luminance":
        gray = 0.299 * rgb_f[..., 0] + 0.587 * rgb_f[..., 1] + 0.114 * rgb_f[..., 2]
    else:
        raise SystemExit(f"unknown he_grayscale.method: {method}")
    return _u8(gray), rgb.astype(np.uint8)


def main() -> int:
    dc = load_config()
    index = load_index(dc)
    slides = index.get("slides", {})
    if not slides:
        log.error("no slides indexed; run 01_index.py first")
        return 1

    target_ds = dc.section("registration").get("thumb_target_downsample", 32)
    ref_cfg = dc.section("reference_channel")
    he_method = dc.section("he_grayscale").get("method", "invert_value")

    for sid, slide in sorted(slides.items()):
        codex = slide.get("codex"); he = slide.get("he")
        if not codex:
            continue
        codex_path = (ROOT / codex["path"]).resolve()
        if not codex_path.exists():
            log.warning("CODEX path missing for %s: %s", sid, codex_path)
            continue
        log.info("CODEX thumb for %s", sid)
        try:
            ref_idx = resolve_reference_channel_index(
                codex.get("channels") or [],
                ref_cfg.get("match_any") or [],
                ref_cfg.get("fallback_index", 0),
            )
            log.info("  ref channel idx=%d (%s)", ref_idx, (codex.get('channels') or [None])[ref_idx])
            if dc.modality in ("codex", "phenocycler", "cycif", "orion"):
                thumb = thumb_codex(codex_path, ref_idx, target_ds)
            else:
                thumb = thumb_via_reader(codex_path, ref_idx, target_ds, dc.modality)
            Image.fromarray(thumb).save(dc.thumbs_dir / f"{sid}__codex.png")
        except Exception as e:
            log.warning("  CODEX thumb failed: %s", e)

        if he:
            he_path = (ROOT / he["path"]).resolve()
            if he_path.exists():
                log.info("H&E thumb for %s", sid)
                try:
                    thumb, rgb = thumb_he(he_path, target_ds, he_method)
                    Image.fromarray(thumb).save(dc.thumbs_dir / f"{sid}__he.png")
                    if rgb is not None:
                        Image.fromarray(rgb).save(dc.thumbs_dir / f"{sid}__he_rgb.png")
                except Exception as e:
                    log.warning("  H&E thumb failed: %s", e)
    log.info("thumbs done → %s", dc.thumbs_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
