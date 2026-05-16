# Tissue H&E + Spatial Proteomics — multiplex preprocessing pipeline

Config-driven preprocessing pipeline for paired tissue H&E + spatial proteomics
datasets across four major modality classes:

| Modality | Reader | Reference channel | Status |
|---|---|---|---|
| **CODEX** (Akoya / Standard BioTools) | `OmeTiffReader` | DAPI | ✅ validated end-to-end |
| **PhenoCycler** | `OmeTiffReader` | DAPI | ✅ same code path as CODEX |
| **CyCIF / t-CyCIF** (Sorger lab) | `OmeTiffReader` | DAPI | ✅ same code path |
| **Orion** (RareCyte CyCIF + same-section H&E) | `OmeTiffReader` + `IdentityStrategy` | DAPI | ✅ identity registration ready |
| **IMC** (Hyperion / Standard BioTools) | `ImcReader` (readimc) | Ir191 / Ir193 (DNA1/DNA2) | ⚠️ implemented to spec, not yet validated against real `.mcd` |
| **MIBI-TOF** (IONpath) | `MibiReader` (per-channel TIFF dir) | dsDNA / Histone H3 | ⚠️ implemented to spec, validated on synthetic ROI |

The reference dataset (`greenwald-gbm-codex`) is a Greenwald 2024 *Cell* GBM
CODEX + paired H&E from the Zenodo 12624860 deposit. New datasets are added by
dropping a YAML config into `configs/datasets/` — no code changes required.

## Pipeline

```
00_unpack    → unpack raw archives (per dataset config)
01_index     → build preproc/<dataset>/index.json (cube + paired H&E metadata)
01b_panels   → resolve channel→marker mapping; verify against OME-XML structure
02_thumbs    → write low-res grayscale thumbs (cube reference channel + H&E 1−V)
03_register  → cross-section LoFTR / ORB / IdentityStrategy → affine
04_patchify  → batched GPU warp + tile to paired (H&E, cube) HDF5
```

Outputs land at `preproc/<dataset>/`. Run with:

```bash
DATASET=greenwald-gbm-codex PY=~/miniforge3/envs/gbm-preproc/bin/python bash run_pipeline.sh
```

To regenerate the preview montage shown in the prompt — an 8 × 8 grid of
`[H&E | DAPI grayscale | H&E + cyan DAPI overlay]` tiles from
`preproc/greenwald-gbm-codex/patches/ZH1041_T1.h5` — run:

```bash
DATASET=greenwald-gbm-codex ~/miniforge3/envs/gbm-preproc/bin/python \
  scripts/render_patch_mosaic.py --sample ZH1041_T1 --channel 0 --max-width 2048
```

The same command is available as:

```bash
make preview PY=~/miniforge3/envs/gbm-preproc/bin/python
```

Standalone preview outputs are written under `results/mosaics/<dataset>/`.
The larger pipeline cache remains under gitignored `preproc/<dataset>/`.

GPU is optional — `04_patchify` falls back to scipy CPU when CUDA isn't
available; `03_register` falls back from kornia LoFTR to skimage ORB. End-to-end
on a single CODEX slide (43 channels, 11k × 11k full-res) is ~80 s on a single
GB10, ~50 min on CPU-only.

## Adding a new dataset

1. Copy a template:
   ```bash
   cp configs/datasets/_template_codex.yaml  configs/datasets/<dataset_id>.yaml
   # or _template_cycif.yaml / _template_imc.yaml / _template_mibi.yaml
   ```
2. Edit `dataset_id`, `paths.raw_subdir`/`preproc_subdir`, the panel source,
   patient-id strip rules, and the reference-channel `match_any` list.
3. Drop raw files under `raw/<dataset_id>/`.
4. Run: `DATASET=<dataset_id> bash run_pipeline.sh`.

The pipeline auto-detects modality from the config and routes to the correct
reader + registration strategy.

## Layout

```
.
├── README.md                                # this file
├── run_pipeline.sh                          # DATASET= env var picks the dataset
├── configs/datasets/                        # production configs — pipeline runs end-to-end
│   ├── greenwald-gbm-codex.yaml             # reference dataset (CODEX + H&E)
│   ├── lin-2022-cycif-crc-3d-atlas.yaml     # 2nd validated config (CyCIF, AWS)
│   ├── lin-2023-orion-crc-sample.yaml       # 3rd validated config (Orion same-section)
│   ├── _template_codex.yaml
│   ├── _template_cycif.yaml
│   ├── _template_imc.yaml
│   └── _template_mibi.yaml
├── configs/_validation/                     # reader smoke-tests — pipeline halts partway by design
│   ├── README.md
│   └── bodenmiller-covid-imc-sample.yaml    # IMC reader validation (no paired H&E in deposit)
├── scripts/
│   ├── _common.py                           # config-aware path constants + helpers
│   ├── dataset_config.py                    # YAML loader + schema
│   ├── readers.py                           # OmeTiffReader / ImcReader / MibiReader
│   ├── registration_strategies.py           # AffineLoFTR / AffineOrb / Identity
│   ├── panels.py                            # named panel registry (Greenwald + future)
│   ├── 00_unpack.py / 01_index.py / 01b_panels.py / 02_thumbs.py
│   ├── 03_register.py / 04_patchify.py
│   ├── mosaics.py / render_patch_mosaic.py  # deterministic H&E/CODEX preview regeneration
│   ├── qc_napari.py                         # interactive single-cell-resolution viewer
│   └── qc_protein_panels.py                 # per-channel patch montages
├── raw/<dataset_id>/                        # downloaded artefacts (gitignored)
└── preproc/<dataset_id>/                    # pipeline outputs (gitignored)
    ├── index.json
    ├── thumbs/
    ├── registration/                        # affine + checkerboard / blend / quiver QC
    ├── patches/<sid>.h5                     # paired patches + verified marker names
    ├── panels/<sid>.csv                     # per-slide channel→marker audit
    └── logs/                                # per-step log files
```

## Reference dataset — Greenwald 2024 *Cell* GBM CODEX

The `greenwald-gbm-codex` dataset is the original target the pipeline was
designed and validated against.

- **Paper:** Greenwald, Galili-Darnell, Hoefflin et al., *Cell* 2024 — "Integrative
  spatial analysis reveals a multi-layered organization of glioblastoma."
- **Deposit:** https://zenodo.org/records/12624860 (CC BY 4.0)
- **Cohort:** 12 GBM CODEX + 6 IDH-mut Visium + 1 GBM organoid CODEX + 13 fresh-frozen GBM Visium
- **claw-data DATA_CARD:** `../claw-data/discovery/multiplex_tissue_imaging/candidates/zenodo-12624860-greenwald-gbm-codex/FINAL.md`
- **Panel:** 43 markers (DAPI + 14 cycles × Atto550/AF488/AF647 with structured gaps in cycles 12-15). Verified against `paper/mmc1.pdf` Table S3 + the deposit's QuPath v2 project per-slide `server.json`.

Running the pipeline against this dataset reproduces:
- `preproc/greenwald-gbm-codex/registration/<sid>.{json,_overlay,_blend,_quiver,_pair}.png`
- `preproc/greenwald-gbm-codex/patches/<sid>.h5` — `(N=256, 512, 512, 3) uint8` H&E + `(N=256, 43, 512, 512) uint8` CODEX, `xy_he` coords, verified marker names in attrs
- `preproc/greenwald-gbm-codex/patches/<sid>_panels/ch00_DAPI.png` … `ch42_CD19.png`

## Validation status per modality

| Modality | Architecture wired | Reader implemented | Real-data validated |
|---|---|---|---|
| CODEX/CyCIF/PhenoCycler/Orion | ✅ | ✅ | ✅ Greenwald end-to-end |
| IMC | ✅ | ✅ (readimc 0.9.2 API) | ⚠️ pending — real `.mcd` not yet downloaded; smoke-tested with synthetic empty file (correct error path) |
| MIBI-TOF | ✅ | ✅ (per-channel TIFF dir) | ⚠️ pending — real ROI not yet downloaded; validated on synthetic 3-channel ROI |

When validating IMC/MIBI: drop a sample under `raw/<dataset_id>/`, write the
config from the template, and `DATASET=<id> bash run_pipeline.sh`. If the
reader hits a real-data quirk not covered by the documented spec, it'll surface
at `02_thumbs` or `04_patchify` with a clean stack trace pointing at the
specific reader method.

## Dependencies

```
tifffile>=2024.7  zarr>=3.0  numpy>=1.26  scipy>=1.11  scikit-image>=0.22
Pillow>=10.0     h5py>=3.10 imagecodecs>=2024.6  pyyaml>=6.0

# GPU (optional — pipeline degrades to CPU when missing)
torch>=2.5  kornia>=0.7

# Modality-specific (optional — only needed for the corresponding modality)
readimc>=0.9   pandas        # IMC (.mcd parsing)
```
