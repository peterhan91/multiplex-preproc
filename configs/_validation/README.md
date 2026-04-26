# `configs/_validation/` — reader smoke-tests, NOT production datasets

Configs in this directory exist to **validate reader / pipeline component
behavior against real data**, not to run the pipeline end-to-end. They
typically lack paired H&E (or some other prerequisite) and the pipeline halts
partway through — by design. The bundled samples here let us catch reader-
binding regressions when changing `scripts/readers.py` or modality-specific
code paths without needing the user to download a full real dataset first.

If you want to **run the full pipeline** on a dataset, use a config from
`configs/datasets/` instead, or copy a `_template_*.yaml` from there.

## Current validation samples

| Config | Purpose | What it exercises | What it skips |
|---|---|---|---|
| `bodenmiller-covid-imc-sample.yaml` | ImcReader against real Hyperion `.mcd` | 00_unpack → 01_index → 01b_panels → 02_thumbs (cube reader, modality dispatch, mcd_metadata panel source, priority-ordered ref channel) | 03_register, 04_patchify (no paired H&E in deposit; IMC is FOV-based not WSI by physics) |

## How to invoke

Same as production datasets:

```bash
DATASET=bodenmiller-covid-imc-sample bash run_pipeline.sh
```

The loader (`scripts/dataset_config.py`) searches both `configs/datasets/`
and `configs/_validation/` so the env-var name doesn't change.

## When to add a new entry here

- You have a small public sample of a modality with no easy paired-H&E source
  (e.g. a single MCD or MIBI ROI directory) and want to keep a known-working
  config so future reader changes can be smoke-tested against real bytes.
- Or: a dataset is in `configs/datasets/` but you want to keep a tiny
  always-present subset for CI-style validation that doesn't need a TB-scale
  download.
