#!/usr/bin/env python3
"""Cross-section / same-section registration of cube ↔ H&E thumbnails.

Strategy is dataset-config driven:

  registration:
    strategy: affine_loftr   # affine_loftr | affine_orb | identity
    thumb_target_downsample: 32
    loftr_max_side: 1024
    ransac_residual_px: 4.0
    min_inliers: 8

The transform maps coordinates in the cube thumbnail space → H&E thumbnail
space; 04_patchify upscales the matrix to full-resolution coordinates.

QC artifacts saved per slide:
  preproc/<dataset>/registration/<sid>.json     transform + diagnostics
  preproc/<dataset>/registration/<sid>_overlay.png   checkerboard
  preproc/<dataset>/registration/<sid>_blend.png     R=H&E, G=warped cube ref
  preproc/<dataset>/registration/<sid>_quiver.png    per-tile residual shifts
  preproc/<dataset>/registration/<sid>_pair.png      side-by-side
  preproc/<dataset>/registration/_summary.json
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
from PIL import Image

from _common import (
    detect_device,
    get_logger,
    load_config,
    load_index,
)
from registration_strategies import make_strategy

log = get_logger("register")
DEV = detect_device()


# ---------------------------------------------------------------- QC

def _hist_match(src: np.ndarray, ref: np.ndarray) -> np.ndarray:
    s = src.ravel(); r = ref.ravel()
    s_vals, s_idx, s_cnt = np.unique(s, return_inverse=True, return_counts=True)
    r_vals, r_cnt = np.unique(r, return_counts=True)
    s_q = np.cumsum(s_cnt).astype(np.float64) / s.size
    r_q = np.cumsum(r_cnt).astype(np.float64) / r.size
    mapped = np.interp(s_q, r_q, r_vals)
    return mapped[s_idx].reshape(src.shape).astype(src.dtype)


def _checkerboard(a: np.ndarray, b: np.ndarray, sz: int = 64) -> np.ndarray:
    h, w = a.shape[:2]
    out = a.copy()
    if out.ndim == 2:
        out = np.stack([out] * 3, axis=-1)
    bb = b
    if bb.ndim == 2:
        bb = np.stack([bb] * 3, axis=-1)
    if bb.shape[:2] != out.shape[:2]:
        from PIL import Image as _Im
        bb = np.asarray(_Im.fromarray(bb).resize((w, h)))
    yy, xx = np.indices((h, w))
    mask = ((yy // sz) + (xx // sz)) % 2 == 0
    out[mask] = bb[mask]
    return out


def _blend_rg(he_gray: np.ndarray, warped_cube: np.ndarray) -> np.ndarray:
    g = _hist_match(warped_cube, he_gray)
    rgb = np.zeros((*he_gray.shape, 3), dtype=np.uint8)
    rgb[..., 0] = he_gray
    rgb[..., 1] = g
    return rgb


def _shift_quiver(he_gray: np.ndarray, warped_cube: np.ndarray, grid: int = 24, win_min: int = 64) -> np.ndarray:
    from skimage.registration import phase_cross_correlation
    h, w = he_gray.shape
    cell_h = max(win_min, h // grid); cell_w = max(win_min, w // grid)
    nrows = h // cell_h; ncols = w // cell_w
    if nrows < 2 or ncols < 2:
        return np.stack([he_gray] * 3, axis=-1)
    shifts = np.zeros((nrows, ncols, 2), dtype=np.float32)
    mags = np.zeros((nrows, ncols), dtype=np.float32)
    for r in range(nrows):
        for c in range(ncols):
            y0, x0 = r * cell_h, c * cell_w
            y1, x1 = y0 + cell_h, x0 + cell_w
            ref = he_gray[y0:y1, x0:x1].astype(np.float32)
            mov = warped_cube[y0:y1, x0:x1].astype(np.float32)
            if ref.std() < 2.0 or mov.std() < 2.0:
                continue
            try:
                shift, _e, _p = phase_cross_correlation(ref, mov, upsample_factor=4, normalization=None)
            except Exception:
                continue
            shifts[r, c] = shift; mags[r, c] = float(np.linalg.norm(shift))

    canvas = np.stack([he_gray] * 3, axis=-1).astype(np.float32)
    if mags.max() > 0:
        mag_img = np.kron(mags, np.ones((cell_h, cell_w), dtype=np.float32))
        pad_h = max(0, h - mag_img.shape[0]); pad_w = max(0, w - mag_img.shape[1])
        if pad_h or pad_w:
            mag_img = np.pad(mag_img, ((0, pad_h), (0, pad_w)))
        mag_img = mag_img[:h, :w]
        m = mag_img / max(1e-6, mag_img.max())
        red = np.zeros_like(canvas); red[..., 0] = 255.0
        canvas = canvas * (1 - 0.45 * m[..., None]) + red * (0.45 * m[..., None])
    out = np.clip(canvas, 0, 255).astype(np.uint8)

    from PIL import Image as _Im, ImageDraw as _Dr
    im = _Im.fromarray(out); dr = _Dr.Draw(im)
    arrow_scale = min(cell_h, cell_w) * 0.6 / max(1.0, mags.max())
    for r in range(nrows):
        for c in range(ncols):
            dy, dx = shifts[r, c]
            if dy == 0 and dx == 0:
                continue
            cx = c * cell_w + cell_w / 2; cy = r * cell_h + cell_h / 2
            ex = cx + dx * arrow_scale; ey = cy + dy * arrow_scale
            dr.line([(cx, cy), (ex, ey)], fill=(255, 255, 0), width=1)
            dr.ellipse([ex - 1, ey - 1, ex + 1, ey + 1], fill=(255, 255, 0))
    return np.asarray(im)


def _warp_skimage(im: np.ndarray, mat: np.ndarray, out_shape) -> np.ndarray:
    from skimage.transform import warp, AffineTransform
    inv = AffineTransform(matrix=np.linalg.inv(mat))
    warped = warp(im.astype(np.float32) / 255.0, inv, output_shape=out_shape, preserve_range=False)
    return (warped * 255).astype(np.uint8)


# ---------------------------------------------------------------- main

def register_one(sid: str, dc) -> dict | None:
    codex_thumb = dc.thumbs_dir / f"{sid}__codex.png"
    he_thumb = dc.thumbs_dir / f"{sid}__he.png"
    if not (codex_thumb.exists() and he_thumb.exists()):
        log.warning("missing thumbnail(s) for %s; skipping", sid)
        return None
    src = np.asarray(Image.open(codex_thumb).convert("L"))
    dst = np.asarray(Image.open(he_thumb).convert("L"))
    log.info("registering %s (cube %s → H&E %s)", sid, src.shape, dst.shape)

    reg_cfg = dc.section("registration")
    strat_name = reg_cfg.get("strategy", "affine_loftr")
    # CPU/GPU dispatch — LoFTR needs CUDA; downgrade to ORB if not available.
    if strat_name == "affine_loftr" and not (DEV.use_cuda and DEV.kornia_available):
        log.info("  no CUDA/kornia → downgrading affine_loftr → affine_orb")
        strat_name = "affine_orb"
    strategy = make_strategy(
        strat_name,
        max_side=reg_cfg.get("loftr_max_side", 1024),
        residual_px=reg_cfg.get("ransac_residual_px", 4.0),
        min_inliers=reg_cfg.get("min_inliers", 8),
    )
    mat, n_matched, n_inliers = strategy.fit(src, dst, device="cuda" if DEV.use_cuda else "cpu")
    if mat is None:
        out = {"sample_id": sid, "ok": False, "matches": n_matched, "inliers": n_inliers,
               "backend": strategy.name}
        (dc.reg_dir / f"{sid}.json").write_text(json.dumps(out, indent=2))
        return out

    warped = _warp_skimage(src, mat, dst.shape)
    Image.fromarray(_checkerboard(dst, warped)).save(dc.reg_dir / f"{sid}_overlay.png")
    Image.fromarray(_blend_rg(dst, warped)).save(dc.reg_dir / f"{sid}_blend.png")
    try:
        Image.fromarray(_shift_quiver(dst, warped)).save(dc.reg_dir / f"{sid}_quiver.png")
    except Exception as e:
        log.warning("  quiver QC failed: %s", e)

    pair = np.zeros((max(src.shape[0], dst.shape[0]), src.shape[1] + dst.shape[1] + 8, 3), dtype=np.uint8)
    pair[: src.shape[0], : src.shape[1]] = np.stack([src] * 3, axis=-1)
    pair[: dst.shape[0], src.shape[1] + 8 :] = np.stack([dst] * 3, axis=-1)
    Image.fromarray(pair).save(dc.reg_dir / f"{sid}_pair.png")

    out = {
        "sample_id": sid, "ok": True,
        "matches": n_matched, "inliers": n_inliers,
        "thumb_shape_codex": list(src.shape),
        "thumb_shape_he": list(dst.shape),
        "affine_thumb_codex_to_he": mat.tolist(),
        "device": str(DEV),
        "backend": strategy.name,
    }
    (dc.reg_dir / f"{sid}.json").write_text(json.dumps(out, indent=2))
    log.info("  saved transform → %s", dc.reg_dir / f"{sid}.json")
    return out


def main() -> int:
    dc = load_config()
    index = load_index(dc)
    slides = index.get("slides", {})
    if not slides:
        log.error("no slides; run 01_index first")
        return 1
    log.info("device: %s  strategy: %s", DEV, dc.section("registration").get("strategy"))
    summary = []
    for sid in sorted(slides):
        slide = slides[sid]
        if not (slide.get("codex") and slide.get("he")):
            log.info("skip %s — missing CODEX or H&E", sid)
            continue
        out = register_one(sid, dc)
        if out:
            summary.append(out)
    (dc.reg_dir / "_summary.json").write_text(json.dumps(summary, indent=2))
    log.info("registration done → %s", dc.reg_dir)
    return 0


if __name__ == "__main__":
    sys.exit(main())
