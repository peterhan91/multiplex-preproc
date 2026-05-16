# Spatial proteomics datasets — catalog

This is the running index of paired-H&E spatial-proteomics datasets we have
either (a) **wired into this repo's pipeline** under `configs/datasets/` and
`configs/_validation/`, or (b) **surveyed for inclusion** via the dataset
discovery system at `~/codes/claw-data/discovery/multiplex_tissue_imaging/`.

The authoritative survey output is the v7 slate (12 candidates that PASSED
all 7 hard filters — modality, panel size, paired H&E, open access, etc.).
The full per-candidate deep cards live at
`claw-data/discovery/multiplex_tissue_imaging/candidates/<slug>/DATA_CARD.md`.
This file mirrors the headlines so the preprocessing-pipeline repo is
self-explanatory.

Survey date: **2026-04-29** (v7, fully-fresh). Re-run discovery in
`claw-data/` to refresh.

---

## Wired into this repo

Datasets with a YAML config and a verified end-to-end run.

| Status | Dataset ID | Modality | Cohort | Panel | Source | Config |
|---|---|---|---|---|---|---|
| ✅ end-to-end validated | `greenwald-gbm-codex` | CODEX (Akoya) | 12 GBM patients (Greenwald 2024 *Cell*) | 43-ch (DAPI + 14 cycles × Atto550/AF488/AF647) | Zenodo [12624860](https://zenodo.org/records/12624860) — CC BY 4.0 | [configs/datasets/greenwald-gbm-codex.yaml](../configs/datasets/greenwald-gbm-codex.yaml) |
| ✅ end-to-end validated | `lin-2022-cycif-crc-3d-atlas` | t-CyCIF | 17 CRC patients (CRC1 3D, CRC2-17 2D) — Lin 2023 *Cell* | 36-ch (HTAN CRC202105) | AWS S3 `lin-2021-crc-atlas` — open anon | [configs/datasets/lin-2022-cycif-crc-3d-atlas.yaml](../configs/datasets/lin-2022-cycif-crc-3d-atlas.yaml) |
| ✅ end-to-end validated | `lin-2023-orion-crc-sample` | Orion (RareCyte CyCIF + same-section H&E) | 1 sample (CRC33_01) — Lin 2023 *Nat Cancer* | 16-18 + Hoechst | AWS S3 `lin-2023-orion-crc/data/CRC33_01/` — open anon | [configs/datasets/lin-2023-orion-crc-sample.yaml](../configs/datasets/lin-2023-orion-crc-sample.yaml) |
| ⚠️ reader-only smoke test (no paired H&E in deposit) | `bodenmiller-covid-imc-sample` | IMC (Hyperion) | 1 MCD × 4 acquisitions × 47 channels — Bodenmiller COVID lung | ActivationPanel 47-ch | Zenodo [4637034](https://zenodo.org/records/4637034) — A19_33_20210121_ActivationPanel.mcd, 632 MB | [configs/_validation/bodenmiller-covid-imc-sample.yaml](../configs/_validation/bodenmiller-covid-imc-sample.yaml) |

How to add another: copy `configs/datasets/_template_<modality>.yaml`,
edit the TODO fields, drop raw files into `raw/<dataset_id>/`, run
`DATASET=<id> bash run_pipeline.sh`.

---

## Surveyed v7 slate — 12 PASS candidates

The 12 datasets that passed all 7 hard filters in the most recent discovery
sweep. Wired-up rows above appear here too with a ✅ in the **Wired** column.
Patient count is conservative (no double-counting overlapping atlas wraps).

| Rank | Wired | Slug | Source | N (cancer/other/healthy/total) | Modality | Panel | H&E location | License |
|---|---|---|---|---|---|---|---|---|
| 1 | — | `mendeley-d87vg86zd8` | Mendeley | 79 / 0 / 0 / **79** | MIBI-TOF | 37-plex | Risom 2022 *Cell* — serial-section + co-registration (PMC8792442) | CC BY 4.0 |
| 2 | — | `lin-2023-orion-crc-cohort1-2` | HTAN / Synapse | 74 / 0 / 0 / **74** | Orion (same-section) | 16-18 + Hoechst | same-section, paper PMC10368530 | NCI Moonshot (open-claimed; file-level watchlist) |
| 3 | — | `hubmap-greenbaum-decidua-2023` | HuBMAP | 0 / 0 / 66 / **66** | MIBI-TOF | 37-plex | Greenbaum 2023 PMC10076224 — implantation site H&E | CC BY 4.0 (conv.) |
| 4 | — | `angelolab-keren-2018-tnbc-mibi` | angelolab.com | 41 / 0 / 0 / **41** | MIBI-TOF | 36-plex | Keren 2018 PMC6132072 — adjacent-section H&E | unstated (research-use) |
| 5 | — | `tcia-crc_ffpe-codex_cellneighs` | TCIA | 35 / 0 / 0 / **35** | CODEX | 56-plex (Schurch 2020) | TCIA collection: "Histopathology, CODEX images" | CC BY 4.0 |
| 6 | ✅ `greenwald-gbm-codex` | `zenodo-12624860` | Zenodo | 19 / 0 / 0 / **19** | CODEX + Visium ST | 43-plex (panel verified in this repo) | `HE_scans.zip` 131 MB in deposit | CC BY 4.0 |
| 7 | — | `figshare-7174914` | Figshare+ | 5 / 0 / 12 / **17** | CODEX (PhenoCycler-Fusion) | 54-plex | sibling articles 25126898 + 25127633 (18 H&E QPTIFFs, same-section) | CC0 |
| 8 | — | `tcia-cptac-glioblastoma-codex` | TCIA | 10 / 0 / 2 / **12** | CODEX (Akoya Keyence) | 25-plex | adjacent-section H&E "after multiplex imaging" verbatim | CC BY 4.0 |
| 9 | — | `wustl-htan-2024-2d3d-ding` | IDC | 10 / 0 / 0 / **≥10** | IMC (DICOM-SM) | ≥35-plex | same-StudyInstanceUID DICOM-SM `H&E` series + `IMC` series | CC BY 4.0 (IDC conv.) |
| 10 | — | `mendeley-dr5fkgtrb6` | Mendeley | 0 / ≥9 / 0 / **≥9** (TB) | MIBI-TOF | 37-plex | McCaffrey 2022 PMC8810384 — paired H&E + MIBI in figures | CC BY 4.0 |
| 11 | — | `hubmap-mibi-uterus-atlas` | HuBMAP | 0 / 0 / ≥66 / **≥66** (atlas wraps Greenbaum + neighbors) | MIBI-TOF | 37-plex | paper PMC10076224 | CC BY 4.0 (conv.) |
| 12 | — | `zenodo-5945388` | Zenodo | 0 / 0 / 0 (archival multi-tissue; 126 ROIs) | MIBI-TOF | 16-plex | NanoZoomer 40× WSI H&E + serial-section co-reg, PMC10357968 | CC BY 4.0 |

**Slate totals (no atlas double-count):** 414 patients across 12 datasets —
233 cancer, 9 other-disease, 80 healthy + 1 method-validation cohort (Risom
archival). Modality split: 8 MIBI-TOF, 3 CODEX (incl. PhenoCycler-Fusion),
1 Orion, 1 IMC. Cancer-only datasets: 7.

Per-candidate deep cards: `~/codes/claw-data/discovery/multiplex_tissue_imaging/candidates/<slug>/DATA_CARD.md`.

---

## Surveyed — excluded at deep-fetch (won't be wired up)

7 candidates that looked viable from search-time metadata but failed a hard
filter after the deep-fetch reading.

| Slug | Failed filter | Why |
|---|---|---|
| `engjen-2025-cycif-tmas-syn50134757` | #7 (paired H&E WSI) | PMC11948582 grep returned zero H&E image-context hits; cycIF_TMAs Synapse deposit is TMA-only |
| `biostudies-S-BIAD1579` | #7 | 4,593-file inventory contains zero H&E / hematoxylin / brightfield / `.svs` / `.ndpi` / `.qptiff` entries |
| `hubmap-codex-atlas` | #7 | Hickey 2023 *Cell* (PMC10356619) only mentions H&E as panel-validation morphology QC — not paired WSI |
| `biostudies-S-BIAD1344` | #4b (panel ≥ 10) | Panel 1 = 4 phenotype markers, Panel 2 = 4 phenotype markers — max single-panel plex = 4 |
| `mendeley-cydmwsfztj` | #7 + N corrected | Damond 2019 paired imaging is IF (not H&E); deposit N=12, not 88 as claimed in search-time metadata |
| `biostudies-S-BIAD2027` | #7 | 144 files all `.mcd` (Hyperion IMC raw) — zero H&E |
| `idr-project-3052` | #7 (deferred) | Greenwald 2026 *Nat Cancer* PMC fulltext not yet available (12-month embargo); IDR has 23 TMAs of MIBI-only |

## Download-broken watchlist (passed filters; bulk download path is awkward)

| Slug | Concern |
|---|---|
| `lin-2023-orion-crc-cohort1-2` | Filter #2 PASS-with-watchlist — HTAN Synapse anonymous probe returned 403; needs JS-rendered probe of HTA7 atlas |
| `hubmap-mibi-uterus-atlas` | Per-asset anonymous HTTPS returns 404; Globus required for bulk transfer |
| `hubmap-greenbaum-decidua-2023` | Same — HuBMAP asset endpoints anonymous-404; needs free Globus account |

---

## How to refresh this catalog

1. Open `~/codes/claw-data/SESSION_STATE.md` — confirms current state of all topics.
2. Re-run discovery on the multiplex_tissue_imaging topic from `claw-data/`
   (the `dx-orchestrator` agent owns this) to produce a v8 slate.
3. Diff the new `DISCOVERY_SUMMARY.md` against this catalog and update both
   sections (wired + surveyed).

Note: claw-data archives prior pipeline lineages as
`discovery/multiplex_tissue_imaging_v<N>_archive_<YYYYMMDD>/` — only the
non-archived `discovery/multiplex_tissue_imaging/` is current.
