# Patch Mosaics

Generated preview mosaics live here, grouped by dataset id. They are derived
from `preproc/<dataset>/patches/*.h5`, so image files in this folder are
ignored by git.

Regenerate the prompt-sized Greenwald CODEX/H&E preview with:

```bash
DATASET=greenwald-gbm-codex ~/miniforge3/envs/gbm-preproc/bin/python \
  scripts/render_patch_mosaic.py --sample ZH1041_T1 --channel 0 --max-width 2048
```
