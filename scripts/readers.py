"""Modality-aware cube readers.

A reader exposes the same interface so 02_thumbs / 04_patchify don't care
about file format:

    reader = open_reader(path, modality, **opts)
    reader.C, reader.H, reader.W            # dimensions
    reader.channel_names: list[str]         # ordered marker labels (when known)
    reader.preload()                        # decompress full cube into RAM
    reader.crop(y, x, h, w) -> (C, H, W) np.ndarray
    reader.close()

Implementations:
  - OmeTiffReader  (CODEX, CyCIF, t-CyCIF, PhenoCycler, Orion — already used)
  - ImcReader      (Hyperion .mcd; per-channel TIFF directories)  — STUB
  - MibiReader     (IONpath per-channel TIFF directories)         — STUB

The stubs raise NotImplementedError with a clear "install X / wait for sample
data" message rather than silently returning empty arrays.
"""
from __future__ import annotations

from pathlib import Path
import logging

import numpy as np
import tifffile

log = logging.getLogger("readers")


# ---------------------------------------------------------------- OME-TIFF

class OmeTiffReader:
    """Reader for CODEX/CyCIF/PhenoCycler-style OME-TIFFs.

    Single multi-channel pyramidal TIFF with axes CYX (or YXC). Channels are
    named in the OME-XML; pyramid level 0 is full resolution.
    """

    modality_class = "fluorescence-cyclic"

    def __init__(self, path: Path, preload: bool = True):
        self.path = Path(path)
        self._tf = tifffile.TiffFile(self.path)
        s0 = self._tf.series[0]
        self.axes = s0.axes
        self.shape = list(s0.shape)
        # CYX is the standard for CODEX/CyCIF cubes; YXC for some platforms;
        # IYX is used by Orion (OME-TIFF where each Channel is stored as a
        # separate Image / IFD). All three layouts present as a 3D ndarray
        # (channels-first or channels-last) and behave identically for crops.
        if self.axes in ("CYX", "IYX"):
            self.C, self.H, self.W = self.shape
            self._axis_order = "CYX"
        elif self.axes == "YXC":
            self.H, self.W, self.C = self.shape
            self._axis_order = "YXC"
        else:
            raise RuntimeError(f"unsupported OME-TIFF axes: {self.axes} (shape={self.shape})")

        # OME-XML channel names
        import re
        ome = self._tf.ome_metadata or ""
        self.channel_names = re.findall(r"<Channel[^>]*Name=\"([^\"]+)\"", ome)

        self._mem = None
        if preload:
            self._do_preload()

    def _do_preload(self):
        import time
        t0 = time.time()
        s0 = self._tf.series[0]
        arr = s0.asarray()
        if self._axis_order == "YXC":
            arr = np.moveaxis(arr, -1, 0)
        self._mem = np.ascontiguousarray(arr)
        log.info(
            "  preloaded OME-TIFF (%s, %.2f GB) in %.1fs",
            self._mem.shape, self._mem.nbytes / 1e9, time.time() - t0,
        )

    def crop(self, y: int, x: int, h: int, w: int) -> np.ndarray:
        y0 = max(0, y); x0 = max(0, x)
        y1 = min(self.H, y + h); x1 = min(self.W, x + w)
        if y1 <= y0 or x1 <= x0:
            return np.zeros((self.C, h, w), dtype=np.float32)
        if self._mem is not None:
            arr = self._mem[:, y0:y1, x0:x1]
        else:
            # Lazy zarr fallback (kept for very large cubes that don't fit RAM)
            import zarr
            z = zarr.open(self._tf.aszarr(level=0), mode="r")
            if self._axis_order == "CYX":
                arr = np.asarray(z[:, y0:y1, x0:x1])
            else:
                arr = np.asarray(z[y0:y1, x0:x1, :]).transpose(2, 0, 1)
        out = np.zeros((self.C, h, w), dtype=arr.dtype)
        out[:, y0 - y : y0 - y + (y1 - y0), x0 - x : x0 - x + (x1 - x0)] = arr
        return out.astype(np.float32)

    def close(self):
        try:
            self._tf.close()
        except Exception:
            pass
        self._mem = None


# ---------------------------------------------------------------- IMC

class ImcReader:
    """Reader for Imaging Mass Cytometry (Hyperion .mcd; readimc-based).

    An MCD file contains a Slide with N Acquisitions (ROIs). Each acquisition
    is its own (C, H, W) image with channel labels = metal isotopes (e.g.
    ``Ir191Di``) and channel names = antibody targets (e.g. ``DNA1``, ``CD45``).

    This reader treats each acquisition as one cube. To pick which acquisition
    pass it via env var ``IMC_ACQUISITION_INDEX`` (default 0); for multi-ROI
    workflows the higher-level pipeline should be extended to iterate.

    Status: implemented against readimc 0.9.2 docs (well-stable API). Not yet
    validated against a real .mcd in this repo — flagged in Phase 3 follow-up.
    """

    modality_class = "mass-cytometry"

    def __init__(self, path: Path, preload: bool = True, acquisition_index: int | None = None):
        try:
            from readimc import MCDFile
        except ImportError as e:
            raise ImportError(
                f"readimc import failed: {e}\n"
                "Install with: pip install readimc pandas\n"
                "(readimc has an undeclared pandas dependency.)\n"
                f"(needed to open {path})"
            ) from e
        import os
        self.path = Path(path)
        if acquisition_index is None:
            acquisition_index = int(os.environ.get("IMC_ACQUISITION_INDEX", "0"))
        self._fh = MCDFile(self.path)
        self._fh.__enter__()
        slides = self._fh.slides
        if not slides:
            raise RuntimeError(f"{self.path}: no slides in MCD")
        # Collect every acquisition from every slide; index across them.
        self._acqs = [acq for s in slides for acq in s.acquisitions]
        if acquisition_index >= len(self._acqs):
            raise IndexError(
                f"{self.path}: acquisition index {acquisition_index} out of range "
                f"(have {len(self._acqs)} acquisitions)"
            )
        acq = self._acqs[acquisition_index]
        self._acq = acq
        # readimc convention (verified 2026-04-26 against a real Bodenmiller
        # COVID MCD): `acq.channel_labels` holds the antibody / target name
        # when stained (e.g. 'HistoneH3', 'DNA1', 'CD45') and a bare metal
        # label (e.g. '80ArAr', '190BCKG') for background channels.
        # `acq.channel_names` holds raw metal isotope codes (e.g. 'In113',
        # 'Ir191'). For pipeline use, prefer the antibody label; fall back to
        # the isotope code only when the label is empty.
        labels = list(acq.channel_labels or [])
        names  = list(acq.channel_names or [])
        self.channel_names = [lbl if lbl else nm for lbl, nm in zip(labels, names)] or names

        if preload:
            arr = self._fh.read_acquisition(acq)  # returns (C, H, W) float32
            self._mem = np.ascontiguousarray(arr).astype(np.float32)
        else:
            self._mem = None
        self.C, self.H, self.W = (self._mem.shape if self._mem is not None
                                  else (len(self.channel_names), int(acq.height_px), int(acq.width_px)))
        self.shape = list((self.C, self.H, self.W))
        self.axes = "CYX"

    def crop(self, y: int, x: int, h: int, w: int) -> np.ndarray:
        if self._mem is None:
            self._mem = np.ascontiguousarray(self._fh.read_acquisition(self._acq)).astype(np.float32)
        y0 = max(0, y); x0 = max(0, x)
        y1 = min(self.H, y + h); x1 = min(self.W, x + w)
        if y1 <= y0 or x1 <= x0:
            return np.zeros((self.C, h, w), dtype=np.float32)
        arr = self._mem[:, y0:y1, x0:x1]
        out = np.zeros((self.C, h, w), dtype=np.float32)
        out[:, y0 - y : y0 - y + (y1 - y0), x0 - x : x0 - x + (x1 - x0)] = arr
        return out

    def close(self):
        try:
            self._fh.__exit__(None, None, None)
        except Exception:
            pass
        self._mem = None


# ---------------------------------------------------------------- MIBI

class MibiReader:
    """Reader for MIBI-TOF data (IONpath per-channel TIFF directory).

    Standard IONpath layout — one ROI per directory, one TIFF per marker:

        <roi_dir>/
            DAPI.tiff       # not always present
            dsDNA.tiff      # nuclear reference (or Histone_H3.tiff)
            CD45.tiff
            Beta-Catenin.tiff
            ...

    All TIFFs in the ROI must share spatial dimensions; channels stack to
    (C, H, W) in alphabetical filename order.

    Status: implemented to the standard IONpath spec but not yet validated
    against a real ROI in this repo — flagged in Phase 3 follow-up.
    """

    modality_class = "mass-cytometry"

    def __init__(self, path: Path, preload: bool = True):
        self.path = Path(path)
        if not self.path.is_dir():
            raise NotADirectoryError(
                f"MibiReader expects a directory of per-channel TIFFs, got: {self.path}"
            )
        files = sorted(self.path.glob("*.tif")) + sorted(self.path.glob("*.tiff"))
        if not files:
            raise RuntimeError(f"{self.path}: no .tif/.tiff files in MIBI ROI directory")
        # de-dup tif vs tiff if both exist (rare)
        seen_stems = set(); uniq = []
        for f in files:
            if f.stem in seen_stems: continue
            seen_stems.add(f.stem); uniq.append(f)
        self._files = uniq
        self.channel_names = [f.stem for f in self._files]
        self.C = len(self._files)

        # Read first file for spatial dims; all channels must match.
        first = tifffile.imread(self._files[0])
        self.H, self.W = first.shape[-2], first.shape[-1]
        self.axes = "CYX"
        self.shape = list((self.C, self.H, self.W))

        if preload:
            arr = np.empty((self.C, self.H, self.W), dtype=np.float32)
            arr[0] = first.astype(np.float32)
            for i, f in enumerate(self._files[1:], 1):
                im = tifffile.imread(f)
                if im.shape[-2:] != (self.H, self.W):
                    raise RuntimeError(
                        f"{f}: shape {im.shape} does not match {self.H}x{self.W}"
                    )
                arr[i] = im.astype(np.float32)
            self._mem = arr
            log.info("  preloaded MIBI ROI (%s, %.2f MB)", self._mem.shape, self._mem.nbytes / 1e6)
        else:
            self._mem = None

    def crop(self, y: int, x: int, h: int, w: int) -> np.ndarray:
        y0 = max(0, y); x0 = max(0, x)
        y1 = min(self.H, y + h); x1 = min(self.W, x + w)
        if y1 <= y0 or x1 <= x0:
            return np.zeros((self.C, h, w), dtype=np.float32)
        if self._mem is not None:
            arr = self._mem[:, y0:y1, x0:x1]
        else:
            arr = np.empty((self.C, y1 - y0, x1 - x0), dtype=np.float32)
            for i, f in enumerate(self._files):
                im = tifffile.imread(f)
                arr[i] = im[y0:y1, x0:x1].astype(np.float32)
        out = np.zeros((self.C, h, w), dtype=np.float32)
        out[:, y0 - y : y0 - y + (y1 - y0), x0 - x : x0 - x + (x1 - x0)] = arr
        return out

    def close(self):
        self._mem = None


# ---------------------------------------------------------------- factory

_BY_MODALITY = {
    "codex":       OmeTiffReader,
    "phenocycler": OmeTiffReader,
    "cycif":       OmeTiffReader,
    "orion":       OmeTiffReader,
    "imc":         ImcReader,
    "mibi":        MibiReader,
}


def open_reader(path: Path, modality: str, preload: bool = True):
    """Return a reader instance for the given file + modality.

    Raises SystemExit if the modality isn't registered (unknown values caught
    by dataset_config validation, but defensive here too).
    """
    cls = _BY_MODALITY.get(modality)
    if cls is None:
        raise SystemExit(f"unknown modality {modality!r} — must be in {sorted(_BY_MODALITY)}")
    return cls(Path(path), preload=preload)


# ---------------------------------------------------------------- helpers

def resolve_reference_channel_index(channel_names: list[str], match_any: list[str], fallback_index: int = 0) -> int:
    """Find the channel whose name matches the highest-priority needle.

    Iterates `match_any` in order and returns the first channel matching the
    first matching needle (substring + case-insensitive). This lets configs
    prefer e.g. ``DNA1`` over ``HistoneH3`` even when both are valid nuclear
    references in the same panel. Falls back to `fallback_index` when nothing
    matches (default 0 = first channel).
    """
    if not channel_names:
        return fallback_index
    for needle in (match_any or []):
        n = needle.lower()
        for i, name in enumerate(channel_names):
            if name and n in name.lower():
                return i
    return fallback_index
