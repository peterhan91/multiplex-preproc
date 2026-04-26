#!/usr/bin/env bash
# Chain the preprocessing pipeline. Each step writes under preproc/<DATASET>/.
#
# Usage:
#   bash run_pipeline.sh                           # default DATASET=greenwald-gbm-codex
#   DATASET=lin-2022-cycif-crc-3d-atlas bash run_pipeline.sh
#   DATASET=greenwald-gbm-codex PY=~/miniforge3/envs/gbm-preproc/bin/python bash run_pipeline.sh
set -euo pipefail
cd "$(dirname "$0")"

DATASET="${DATASET:-greenwald-gbm-codex}"
PY="${PY:-python3}"

echo "▸ dataset: $DATASET"
echo "▸ python : $PY ($($PY -V 2>&1))"

export DATASET

echo "▸ 00_unpack"
$PY scripts/00_unpack.py

echo "▸ 01_index"
$PY scripts/01_index.py

echo "▸ 01b_panels"
$PY scripts/01b_panels.py

echo "▸ 02_thumbs"
$PY scripts/02_thumbs.py

echo "▸ 03_register"
$PY scripts/03_register.py

echo "▸ 04_patchify"
$PY scripts/04_patchify.py

echo "✅ pipeline done — outputs at preproc/$DATASET/"
