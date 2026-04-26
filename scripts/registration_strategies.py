"""Registration strategies — pluggable per dataset.

Three concrete strategies:
  - AffineLoFTRStrategy   : cross-section thumbnail-level LoFTR/skimage ORB → 6-DoF affine
  - IdentityStrategy      : same-section data (Orion, restained mIF+H&E) — affine = I
  - AffineOrbStrategy     : skimage ORB only (CPU baseline; for boxes without GPU)

Each strategy returns a 3×3 homogeneous matrix mapping CODEX-thumb → H&E-thumb
coordinates plus diagnostic counts. 03_register dispatches on the config.
"""
from __future__ import annotations

import logging

import numpy as np

log = logging.getLogger("registration_strategies")


# ---------------------------------------------------------------- shared

def _ransac_affine(src_xy: np.ndarray, dst_xy: np.ndarray, residual_px: float, min_inliers: int):
    from skimage.measure import ransac
    from skimage.transform import AffineTransform
    if len(src_xy) < 12:
        return None, len(src_xy), 0
    model, inliers = ransac(
        (src_xy, dst_xy), AffineTransform,
        min_samples=3, residual_threshold=residual_px, max_trials=2000,
    )
    n = int(inliers.sum()) if inliers is not None else 0
    log.info("  RANSAC inliers: %d / %d", n, len(src_xy))
    if model is None or n < min_inliers:
        return None, len(src_xy), n
    return model.params, len(src_xy), n


# ---------------------------------------------------------------- LoFTR

class AffineLoFTRStrategy:
    name = "affine_loftr"

    def __init__(self, max_side: int = 1024, residual_px: float = 4.0, min_inliers: int = 8):
        self.max_side = max_side; self.residual_px = residual_px; self.min_inliers = min_inliers

    def fit(self, src: np.ndarray, dst: np.ndarray, device: str = "cuda"):
        import torch
        from kornia.feature import LoFTR
        from PIL import Image as _Im

        def _prep(im, scale):
            h, w = im.shape
            if scale < 1.0:
                nh, nw = int(round(h * scale)), int(round(w * scale))
                im = np.asarray(_Im.fromarray(im).resize((nw, nh), _Im.BILINEAR))
            nh8, nw8 = (im.shape[0] // 8) * 8, (im.shape[1] // 8) * 8
            im = im[:nh8, :nw8]
            t = torch.from_numpy(im.astype(np.float32) / 255.0).to(device)
            return t.unsqueeze(0).unsqueeze(0), 1.0 / scale

        s = min(1.0, self.max_side / max(src.shape + dst.shape))
        src_t, src_back = _prep(src, s); dst_t, dst_back = _prep(dst, s)
        matcher = LoFTR(pretrained="outdoor").to(device).eval()
        with torch.no_grad():
            out = matcher({"image0": src_t, "image1": dst_t})
        mkpts0 = out["keypoints0"].cpu().numpy()
        mkpts1 = out["keypoints1"].cpu().numpy()
        conf = out["confidence"].cpu().numpy()
        log.info("  matched %d LoFTR pairs (mean conf %.3f)", len(mkpts0), float(conf.mean()) if len(conf) else 0.0)
        src_xy = mkpts0.astype(np.float64) * src_back
        dst_xy = mkpts1.astype(np.float64) * dst_back
        return _ransac_affine(src_xy, dst_xy, self.residual_px, self.min_inliers)


# ---------------------------------------------------------------- ORB CPU

class AffineOrbStrategy:
    name = "affine_orb"

    def __init__(self, n_kp: int = 2000, residual_px: float = 4.0, min_inliers: int = 8):
        self.n_kp = n_kp; self.residual_px = residual_px; self.min_inliers = min_inliers

    def fit(self, src: np.ndarray, dst: np.ndarray, device: str = "cpu"):
        from skimage.feature import ORB, match_descriptors
        orb = ORB(n_keypoints=self.n_kp, fast_threshold=0.05)
        orb.detect_and_extract(src.astype(np.float32) / 255.0)
        kp_src, des_src = orb.keypoints, orb.descriptors
        orb.detect_and_extract(dst.astype(np.float32) / 255.0)
        kp_dst, des_dst = orb.keypoints, orb.descriptors
        matches = match_descriptors(des_src, des_dst, cross_check=True)
        log.info("  matched %d ORB pairs", len(matches))
        if len(matches) < 12:
            return None, len(matches), 0
        sxy = kp_src[matches[:, 0]][:, ::-1]; dxy = kp_dst[matches[:, 1]][:, ::-1]
        return _ransac_affine(sxy, dxy, self.residual_px, self.min_inliers)


# ---------------------------------------------------------------- Identity

class IdentityStrategy:
    """For datasets where CODEX/MIF and H&E come from the SAME tissue section
    (Orion CyCIF, restained mIF+H&E). Registration is the identity transform
    plus a uniform scale that maps the CODEX-thumb pixel grid to the H&E-thumb
    pixel grid (the two thumbnails may have different shapes when the underlying
    full-res images differ in resolution)."""
    name = "identity"

    def __init__(self, **_ignored):
        pass

    def fit(self, src: np.ndarray, dst: np.ndarray, device: str = "cpu"):
        h_s, w_s = src.shape; h_d, w_d = dst.shape
        # diagonal scale CODEX-thumb → H&E-thumb
        sx = w_d / max(1, w_s); sy = h_d / max(1, h_s)
        mat = np.array([[sx, 0.0, 0.0],
                        [0.0, sy, 0.0],
                        [0.0, 0.0, 1.0]], dtype=np.float64)
        log.info("  identity strategy (same-section): scale x=%.4f y=%.4f", sx, sy)
        return mat, 0, 0


# ---------------------------------------------------------------- factory

_BY_NAME = {
    "affine_loftr": AffineLoFTRStrategy,
    "affine_orb":   AffineOrbStrategy,
    "identity":     IdentityStrategy,
}


def make_strategy(name: str, **opts):
    cls = _BY_NAME.get(name)
    if cls is None:
        raise SystemExit(f"unknown registration strategy {name!r} — must be in {sorted(_BY_NAME)}")
    return cls(**opts)
