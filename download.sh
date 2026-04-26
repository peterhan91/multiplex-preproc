#!/usr/bin/env bash
# Idempotent, resumable fetcher for the Greenwald 2024 Cell GBM CODEX subset.
# Source: Zenodo record 12624860 (CC BY 4.0)
# Selected files total ~3.8 GB.

set -euo pipefail
cd "$(dirname "$0")"
mkdir -p raw

ZENODO_RECORD=12624860
BASE="https://zenodo.org/api/records/${ZENODO_RECORD}/files"

# file_name : declared_size_bytes  (sourced from the cached Zenodo file list)
FILES=(
  "HE_scans.zip:131169118"
  "ZH1041_T1_B.ome.tif:1883747801"
  "CODEX_organoid_qupath_project.zip:1763800621"
  "GBM_ZH1019T1.tar.gz:23447164"
)

total_expected=0
for entry in "${FILES[@]}"; do
  total_expected=$((total_expected + ${entry##*:}))
done
echo "▸ Subset = ${#FILES[@]} files, ${total_expected} bytes (~$((total_expected / 1024 / 1024)) MiB)"

for entry in "${FILES[@]}"; do
  name="${entry%%:*}"
  expected="${entry##*:}"
  out="raw/${name}"

  if [[ -f "$out" ]] && [[ "$(stat -c%s "$out")" == "$expected" ]]; then
    echo "  ✓ ${name} (already present, ${expected} bytes)"
    continue
  fi

  echo "▸ Downloading ${name} (~$((expected / 1024 / 1024)) MiB)"
  curl --location --fail --continue-at - \
       --output "$out" \
       "${BASE}/${name}/content"

  actual="$(stat -c%s "$out")"
  if [[ "$actual" != "$expected" ]]; then
    echo "  ✗ size mismatch: got ${actual}, expected ${expected}" >&2
    exit 1
  fi
  echo "  ✓ ${name} (${actual} bytes)"
done

echo
echo "✅ Subset download complete. Files under raw/:"
ls -lh raw/
