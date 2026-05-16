# Claude project notes — multiplex-preproc

Orientation for any AI assistant editing this repo. Human contributors should
read `README.md` first; this file documents conventions an AI is likely to
need that aren't obvious from the code.

## What this repo is

Config-driven preprocessing pipeline for paired tissue H&E + spatial
proteomics (CODEX / PhenoCycler / CyCIF / Orion / IMC / MIBI-TOF). Goal:
emit HDF5 patches of co-registered `[H&E_RGB, multiplex_cube]` for downstream
ML at full WSI resolution.

## Data flow (one slide)

```
raw/<dataset>/<cube.ome.tif|.mcd|...>            # downloaded once
raw/<dataset>/<paired_he.ome.tif|...>
        ↓ 00_unpack    (extracts QuPath, Visium, HE archives per config.unpack)
preproc/<dataset>/he_raw/                         # extracted H&E + any sidecars
        ↓ 01_index     → preproc/<dataset>/index.json
        ↓ 01b_panels   → resolves channel→marker (qupath_v2 / csv / inline / constant / mcd_metadata)
        ↓ 02_thumbs    → preproc/<dataset>/thumbs/<sid>__{codex,he,he_rgb}.png
        ↓ 03_register  → preproc/<dataset>/registration/<sid>.{json,_overlay.png,_blend.png,_quiver.png,_pair.png}
        ↓ 04_patchify  → preproc/<dataset>/patches/<sid>.h5  (+ _preview.png triptych)
results/mosaics/<dataset>/<sid>_ch<NN>_preview.png   # standalone via render_patch_mosaic.py
```

Every stage is independently re-runnable; downstream stages read the index +
the upstream output PNG/JSON.

## Where each modality is handled

Adding/removing a modality touches **all four** of these (use grep, don't
trust memory):

| File | What changes per modality |
|---|---|
| `scripts/readers.py` | `_BY_MODALITY` factory + reader class (`OmeTiffReader` for CODEX/PhenoCycler/CyCIF/Orion; `ImcReader` for IMC; `MibiReader` for MIBI) |
| `scripts/registration_strategies.py` | `_BY_NAME` registry (`affine_loftr` cross-section, `affine_orb` CPU fallback, `identity` same-section/Orion) |
| `scripts/dataset_config.py::_validate` | `valid_modalities` set |
| `scripts/01_index.py`, `02_thumbs.py`, `04_patchify.py` | OME-TIFF vs `open_reader(...)` branching on `dc.modality` |

Per-modality oddities:
- **CODEX / PhenoCycler / CyCIF**: all `OmeTiffReader` (single multi-channel pyramidal OME-TIFF, CYX/YXC).
- **Orion**: same-section H&E. Use `registration.strategy: identity` — NEVER LoFTR/ORB. Cube and H&E are pixel-aligned by acquisition.
- **IMC**: `.mcd` file, multiple acquisitions per slide. Pick one via `IMC_ACQUISITION_INDEX` env var (default 0). Panel comes from MCD metadata (`panel.source: mcd_metadata`). Channel naming convention: `acq.channel_labels` = antibody target (preferred), `acq.channel_names` = metal isotope (fallback).
- **MIBI**: per-channel TIFF directory, one TIFF per marker. Reader stacks alphabetically. Nuclear ref is `dsDNA` / `Histone_H3` (NOT DAPI).

## Config conventions (gotchas)

- **Active dataset** selected via `DATASET` env var. Defaults to `greenwald-gbm-codex` (`DEFAULT_DATASET` in `scripts/dataset_config.py`).
- **CSV panel paths** support `raw://<rel>`, `preproc://<rel>`, `/abs/path`, or bare-relative (legacy → resolved relative to `preproc/<dataset>/`).
- **H&E discovery** matches sample_id by stem-substring (underscores stripped, case-insensitive). For datasets where cube + H&E share an id prefix (Lin 2022 Orion: `CRC02.ome.tif` + `CRC02-HE.ome.tif`), use `index.he_filename_re` to disambiguate.
- **`he_root_base`**: `preproc` (default) means H&E was extracted from an archive into `preproc/<dataset>/he_raw/`. `raw` means H&E lives in `raw/<dataset>/` (typical for direct S3 pulls like Lin Orion / Bodenmiller).
- **Validation-only configs** live in `configs/_validation/`. The loader searches both `configs/datasets/` and `configs/_validation/` for the same env-var name. See `configs/_validation/README.md` for the contract.

## Reference-channel resolution

`scripts/readers.py::resolve_reference_channel_index` iterates
`config.reference_channel.match_any` in order — **the list is priority-ordered**.
For IMC put `DNA1` before `HistoneH3` even though both are valid, so panels
where both exist resolve consistently.

## GPU vs. CPU dispatch

- LoFTR registration: GPU-required (kornia). Falls back to skimage ORB on CPU automatically (`03_register.py` line ~155).
- Patchify warp: GPU-accelerated `F.grid_sample`, falls back to per-channel scipy `affine_transform` on CPU.
- Probe: `scripts/_common.py::detect_device()` — guards `import torch` / `import kornia` in try/except so the import is safe on CPU-only machines.

## QA / preview rendering

- Per-slide preview written by `04_patchify` itself: `preproc/<dataset>/patches/<sid>_preview.png` (8×8 triptych grid).
- Standalone CLI: `scripts/render_patch_mosaic.py` — same triptych, configurable channel/cols/overlay color. Output goes to `results/mosaics/<dataset>/`.
- Both call into `scripts/mosaics.py` — keep that as the single source of truth for the H&E/channel/overlay rendering. Don't duplicate the to_u8 / overlay logic anywhere else.

## Things that LOOK refactor-bait but aren't

- `_common.py` runs `load_config()` + `ensure_dirs()` at module-import time so the output dirs exist before any pipeline stage touches them. Don't move this — many scripts rely on the side effect.
- Module-level path constants (`RAW`, `PREPROC`, `HE_RAW`, `THUMBS`, etc.) in `_common.py` are back-compat shims for older scripts. New code should use `load_config().<path>` directly; don't delete the shims unless you've grepped every script.

## Where the datasets come from

`docs/datasets.md` — index of what's wired up + the broader survey output
from `~/codes/claw-data/discovery/multiplex_tissue_imaging/`. Refresh by
re-running the discovery pipeline in `claw-data/`, not by hand-editing.

## Conventions for AI edits

- One commit per logical change. Don't bundle unrelated refactors.
- Don't add backwards-compat shims for code you delete in the same PR.
- Don't add tests just because a function looks testable — match existing testing posture (currently: smoke-test imports in CI, no unit tests).
- Don't write comments restating what the code does. Only document non-obvious WHY.
- If you're adding a new modality, the four files in the "modality handled" table above ALL need edits; missing one breaks `01_index.py`.
