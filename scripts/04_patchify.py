#!/usr/bin/env python3
"""Generate paired (H&E, multiplex-cube) patches at full resolution.

Strategy
--------
1. Read the per-sample affine transform from preproc/<dataset>/registration/<sid>.json
   (operates in cube-thumb → H&E-thumb coordinates).
2. Build the full-res analogue by composing pyramid-level scale matrices.
3. Tile the H&E full-res grid into PATCH_PX × PATCH_PX patches (after a tissue
   mask filter — Otsu on H&E grayscale + min-tissue-fraction).
4. For each H&E patch, compute the corresponding cube FOV by inverse-warping
   the patch corners through the affine, then crop a cube patch and resample
   onto the H&E pixel grid.
   - GPU path (when CUDA + torch present): batched F.grid_sample on the entire
     channel cube, chunked across patches.
   - CPU fallback: scipy.ndimage.affine_transform per channel.
5. Persist as HDF5: preproc/<dataset>/patches/<sid>.h5

The cube reader is selected by the dataset's modality (CODEX/CyCIF/IMC/MIBI)
via scripts/readers.open_reader.
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import h5py
import numpy as np
from PIL import Image

from _common import (
    ROOT,
    detect_device,
    get_logger,
    load_config,
    load_index,
)
from readers import open_reader

log = get_logger("patchify")
DEV = detect_device()

PREVIEW_GRID = 8


# ----------------------------------------------------------------- transforms

def _scale_matrix(s: float) -> np.ndarray:
    return np.array([[s, 0, 0], [0, s, 0], [0, 0, 1]], dtype=np.float64)


def upscale_affine(mat_thumb: np.ndarray, ds_codex: float, ds_he: float) -> np.ndarray:
    return _scale_matrix(ds_he) @ mat_thumb @ _scale_matrix(1.0 / ds_codex)


def _codex_bbox(he_xy, patch_px, inv):
    x_he, y_he = he_xy
    he_corners = np.array(
        [[x_he, y_he, 1], [x_he + patch_px, y_he, 1],
         [x_he + patch_px, y_he + patch_px, 1], [x_he, y_he + patch_px, 1]],
        dtype=np.float64,
    ).T
    codex_corners = (inv @ he_corners)[:2].T
    x0, y0 = np.floor(codex_corners.min(axis=0)).astype(int) - 4
    x1, y1 = np.ceil(codex_corners.max(axis=0)).astype(int) + 4
    return int(x0), int(y0), int(x1), int(y1)


# ----------------------------------------------------------------- H&E reader

class HeWsi:
    def __init__(self, path: Path, mode: str = "preload"):
        import tifffile
        self.path = path
        self._tf = tifffile.TiffFile(self.path)
        s0 = self._tf.series[0]
        self.axes = s0.axes
        self.shape = list(s0.shape)
        if self.axes.endswith(("YXC", "YXS")):
            self.H, self.W, self.Cc = self.shape[-3:]; self._mode = "YXC"
        elif self.axes.endswith("CYX"):
            self.Cc, self.H, self.W = self.shape[-3:]; self._mode = "CYX"
        elif self.axes in {"YX"}:
            self.H, self.W = self.shape; self.Cc = 1; self._mode = "YX"
        else:
            raise RuntimeError(f"unsupported H&E axes: {self.axes} (shape={self.shape})")

        self.preload = mode == "preload"
        if self.preload:
            t0 = time.time()
            arr = s0.asarray()
            if self._mode == "CYX":
                arr = np.moveaxis(arr, -3, -1)[..., :3]
            elif self._mode == "YXC":
                arr = arr[..., :3]
            else:
                arr = np.stack([arr] * 3, axis=-1)
            if arr.dtype != np.uint8:
                m = max(1, arr.max())
                arr = (arr / m * 255).astype(np.uint8) if m > 0 else arr.astype(np.uint8)
            self._mem = np.ascontiguousarray(arr)
            log.info("  preloaded H&E (%s, %.2f GB) in %.1fs",
                     self._mem.shape, self._mem.nbytes / 1e9, time.time() - t0)
            self._z = None
        else:
            import zarr
            self._mem = None
            self._z = zarr.open(self._tf.aszarr(level=0), mode="r")

    def crop(self, y, x, h, w) -> np.ndarray:
        y0 = max(0, y); x0 = max(0, x)
        y1 = min(self.H, y + h); x1 = min(self.W, x + w)
        if y1 <= y0 or x1 <= x0:
            return np.zeros((h, w, 3), dtype=np.uint8)
        if self._mem is not None:
            arr = self._mem[y0:y1, x0:x1]
        elif self._mode == "YXC":
            arr = np.asarray(self._z[y0:y1, x0:x1, :3])
        elif self._mode == "CYX":
            arr = np.asarray(self._z[:3, y0:y1, x0:x1]).transpose(1, 2, 0)
        else:
            g = np.asarray(self._z[y0:y1, x0:x1])
            arr = np.stack([g, g, g], axis=-1)
        out = np.zeros((h, w, 3), dtype=arr.dtype)
        out[y0 - y : y0 - y + (y1 - y0), x0 - x : x0 - x + (x1 - x0), :] = arr
        if out.dtype != np.uint8:
            out = (out / max(1, out.max()) * 255).astype(np.uint8) if out.max() > 0 else out.astype(np.uint8)
        return out

    def close(self):
        try:
            self._tf.close()
        except Exception:
            pass
        self._mem = None


# ----------------------------------------------------------------- cube warps

def warp_codex_patch_cpu(cube, he_xy, patch_px, affine_full):
    from scipy.ndimage import affine_transform
    inv = np.linalg.inv(affine_full)
    x_he, y_he = he_xy
    x0, y0, x1, y1 = _codex_bbox(he_xy, patch_px, inv)
    crop = cube.crop(y0, x0, y1 - y0, x1 - x0)
    M = inv[:2, :2][[1, 0]][:, [1, 0]]
    offset_full = inv[:2, 2][[1, 0]]
    offset = M @ np.array([y_he, x_he]) + offset_full - np.array([y0, x0])
    out = np.empty((crop.shape[0], patch_px, patch_px), dtype=np.float32)
    for c in range(crop.shape[0]):
        out[c] = affine_transform(crop[c], matrix=M, offset=offset,
                                  output_shape=(patch_px, patch_px),
                                  order=1, mode="constant", cval=0.0)
    return out


def warp_codex_batch_gpu(cube, xy_list, patch_px, affine_full, device="cuda", chunk=32):
    import torch
    import torch.nn.functional as F
    inv = np.linalg.inv(affine_full)
    inv_t = torch.from_numpy(inv).to(device=device, dtype=torch.float32)
    P = patch_px
    N = len(xy_list)
    out = np.empty((N, cube.C, P, P), dtype=np.float32)
    yo = torch.arange(P, device=device, dtype=torch.float32)
    xo = torch.arange(P, device=device, dtype=torch.float32)
    yo_g, xo_g = torch.meshgrid(yo, xo, indexing="ij")
    with torch.no_grad():
        for s in range(0, N, chunk):
            sub = xy_list[s : s + chunk]
            B = len(sub)
            bboxes = [_codex_bbox(xy, P, inv) for xy in sub]
            crops_np = [cube.crop(y0, x0, y1 - y0, x1 - x0) for (x0, y0, x1, y1) in bboxes]
            h_max = max(c.shape[1] for c in crops_np)
            w_max = max(c.shape[2] for c in crops_np)
            stacked = np.zeros((B, cube.C, h_max, w_max), dtype=np.float32)
            for i, c in enumerate(crops_np):
                stacked[i, :, : c.shape[1], : c.shape[2]] = c
            src = torch.from_numpy(stacked).to(device=device, dtype=torch.float32)
            grids = torch.empty(B, P, P, 2, device=device, dtype=torch.float32)
            for i, (xy, (x0, y0, x1, y1)) in enumerate(zip(sub, bboxes)):
                x_he, y_he = xy
                he_full_x = xo_g + float(x_he); he_full_y = yo_g + float(y_he)
                co_full_x = inv_t[0, 0] * he_full_x + inv_t[0, 1] * he_full_y + inv_t[0, 2]
                co_full_y = inv_t[1, 0] * he_full_x + inv_t[1, 1] * he_full_y + inv_t[1, 2]
                co_crop_x = co_full_x - float(x0); co_crop_y = co_full_y - float(y0)
                grids[i, :, :, 0] = 2.0 * co_crop_x / max(1, w_max - 1) - 1.0
                grids[i, :, :, 1] = 2.0 * co_crop_y / max(1, h_max - 1) - 1.0
            warped = F.grid_sample(src, grids, mode="bilinear",
                                   padding_mode="zeros", align_corners=True)
            out[s : s + B] = warped.cpu().numpy()
    return out


# ----------------------------------------------------------------- tissue mask

def tissue_mask_from_thumb(he_thumb_path: Path) -> np.ndarray:
    from skimage.filters import threshold_otsu
    img = np.asarray(Image.open(he_thumb_path).convert("L"))
    t = threshold_otsu(img)
    return img > t


# ----------------------------------------------------------------- main

def _ds_factor(full_shape, thumb_shape) -> float:
    return full_shape[-2] / thumb_shape[0]


def patchify_one(sid: str, slide_meta: dict, dc) -> dict:
    reg_path = dc.reg_dir / f"{sid}.json"
    if not reg_path.exists():
        return {"sample_id": sid, "ok": False, "reason": "no_registration"}
    reg = json.loads(reg_path.read_text())
    if not reg.get("ok"):
        return {"sample_id": sid, "ok": False, "reason": "registration_failed"}

    codex_path = (ROOT / slide_meta["codex"]["path"]).resolve()
    he_path = (ROOT / slide_meta["he"]["path"]).resolve()
    he_thumb = dc.thumbs_dir / f"{sid}__he.png"

    cube = open_reader(codex_path, dc.modality, preload=True)
    he = HeWsi(he_path)

    pcfg = dc.section("patchify")
    PATCH_PX = pcfg.get("patch_px", 512)
    STRIDE_PX = pcfg.get("stride_px", 512)
    TISSUE_FRAC_MIN = pcfg.get("tissue_frac_min", 0.35)
    MAX_PATCHES = pcfg.get("max_patches_per_slide", 256)
    storage = pcfg.get("storage", {})
    he_dtype = storage.get("he_dtype", "uint8")
    codex_dtype = storage.get("codex_dtype", "uint8")
    compression = storage.get("compression", "lzf")

    ds_codex = _ds_factor((cube.H, cube.W), reg["thumb_shape_codex"])
    ds_he = _ds_factor((he.H, he.W), reg["thumb_shape_he"])
    affine_full = upscale_affine(np.array(reg["affine_thumb_codex_to_he"]), ds_codex, ds_he)
    log.info("  full-res affine ds_codex=%.3f ds_he=%.3f", ds_codex, ds_he)

    mask = tissue_mask_from_thumb(he_thumb)
    mask_h, mask_w = mask.shape
    he_h, he_w = he.H, he.W
    sx, sy = he_w / mask_w, he_h / mask_h

    coords = []
    for y in range(0, he_h - PATCH_PX + 1, STRIDE_PX):
        for x in range(0, he_w - PATCH_PX + 1, STRIDE_PX):
            mx0, my0 = int(x / sx), int(y / sy)
            mx1, my1 = int((x + PATCH_PX) / sx), int((y + PATCH_PX) / sy)
            if my1 <= my0 or mx1 <= mx0:
                continue
            tissue = mask[my0:my1, mx0:mx1].mean()
            if tissue >= TISSUE_FRAC_MIN:
                coords.append((x, y, float(tissue)))
    coords.sort(key=lambda t: -t[2])
    coords = coords[:MAX_PATCHES]
    log.info("  %d candidate patches after tissue filter", len(coords))
    if not coords:
        cube.close(); he.close()
        return {"sample_id": sid, "ok": False, "reason": "no_tissue_patches"}

    out_path = dc.patches_dir / f"{sid}.h5"
    n = len(coords); p = PATCH_PX
    channel_names = (
        slide_meta["codex"].get("markers")
        or slide_meta["codex"].get("channels")
        or [f"ch{i}" for i in range(cube.C)]
    )
    if len(channel_names) != cube.C:
        channel_names = [f"ch{i}" for i in range(cube.C)]

    with h5py.File(out_path, "w") as f:
        g = f.create_group("patches")
        d_he = g.create_dataset("he", (n, p, p, 3), dtype=he_dtype,
                                chunks=(1, p, p, 3), compression=compression)
        d_co = g.create_dataset("codex", (n, cube.C, p, p), dtype=codex_dtype,
                                chunks=(1, cube.C, p, p), compression=compression)
        d_xy = g.create_dataset("xy_he", (n, 2), dtype="int32")
        f.attrs["sample_id"] = sid
        f.attrs["dataset_id"] = dc.dataset_id
        f.attrs["modality"] = dc.modality
        f.attrs["patch_px"] = p
        f.attrs["stride_px"] = STRIDE_PX
        f.attrs["affine_full"] = affine_full
        f.attrs["channels"] = np.array(channel_names, dtype=h5py.string_dtype())

        previews = []
        t0 = time.time()
        if DEV.use_cuda and DEV.torch_available:
            xy_only = [(x, y) for (x, y, _) in coords]
            log.info("    batched GPU warp: %d patches", n)
            co_all_f32 = warp_codex_batch_gpu(cube, xy_only, p, affine_full, device="cuda", chunk=32)
            t_warp = time.time() - t0
            log.info("    GPU warp done in %.1fs (%.1f patches/s)", t_warp, n / max(1e-3, t_warp))
            co_all = (np.clip(co_all_f32 + 0.5, 0, 255).astype(np.uint8)
                      if codex_dtype == "uint8" else co_all_f32.astype(codex_dtype))
            del co_all_f32
            t_w0 = time.time()
            for i, (x, y, _) in enumerate(coords):
                he_patch = he.crop(y, x, p, p)
                d_he[i] = he_patch; d_co[i] = co_all[i]; d_xy[i] = (x, y)
                if i < PREVIEW_GRID * PREVIEW_GRID:
                    previews.append((he_patch, co_all[i]))
                if i % 64 == 0 and i > 0:
                    rate = i / max(1e-3, time.time() - t_w0)
                    log.info("    write %d/%d (%.1f patches/s)", i, n, rate)
        else:
            for i, (x, y, _) in enumerate(coords):
                he_patch = he.crop(y, x, p, p)
                co_patch = warp_codex_patch_cpu(cube, (x, y), p, affine_full)
                if codex_dtype == "uint8":
                    co_patch = np.clip(co_patch + 0.5, 0, 255).astype(np.uint8)
                d_he[i] = he_patch; d_co[i] = co_patch; d_xy[i] = (x, y)
                if i < PREVIEW_GRID * PREVIEW_GRID:
                    previews.append((he_patch, co_patch))
                if i % 32 == 0:
                    rate = (i + 1) / max(1e-3, (time.time() - t0))
                    log.info("    patch %d/%d (%.1f patches/s)", i, n, rate)

    # 3-up preview (H&E | nuclear-ref | overlay) — use channel 0 as ref proxy
    def _to_u8(a):
        a = a.astype(np.float32)
        lo, hi = np.percentile(a, [1, 99.5])
        if hi <= lo: hi = lo + 1.0
        return (np.clip((a - lo) / (hi - lo), 0, 1) * 255).astype(np.uint8)
    rows = []
    side = max(1, int(np.sqrt(min(len(previews), PREVIEW_GRID * PREVIEW_GRID))))
    for r in range(side):
        row = []
        for c in range(side):
            idx = r * side + c
            if idx >= len(previews): break
            he_p, co_p = previews[idx]
            d_u8 = _to_u8(co_p[0])
            d_rgb = np.stack([np.zeros_like(d_u8), d_u8, d_u8], axis=-1)  # cyan
            blend = np.clip(he_p.astype(np.uint16) + (d_rgb.astype(np.uint16) * 0.6).astype(np.uint16), 0, 255).astype(np.uint8)
            row.append(np.hstack([he_p, np.stack([d_u8] * 3, axis=-1), blend]))
        if row: rows.append(np.hstack(row))
    if rows:
        Image.fromarray(np.vstack(rows)).save(dc.patches_dir / f"{sid}_preview.png")

    cube.close(); he.close()
    return {"sample_id": sid, "ok": True, "n_patches": n,
            "h5": str(out_path), "channels": channel_names}


def main() -> int:
    dc = load_config()
    index = load_index(dc)
    slides = index.get("slides", {})
    summary = []
    for sid in sorted(slides):
        slide = slides[sid]
        if not (slide.get("codex") and slide.get("he")):
            continue
        log.info("patchify %s", sid)
        try:
            out = patchify_one(sid, slide, dc)
        except Exception as e:
            log.exception("patchify failed for %s: %s", sid, e)
            out = {"sample_id": sid, "ok": False, "reason": f"exception: {e}"}
        summary.append(out)
    (dc.patches_dir / "_summary.json").write_text(json.dumps(summary, indent=2))
    log.info("patchify done — %d sample(s)", len(summary))
    return 0


if __name__ == "__main__":
    sys.exit(main())
