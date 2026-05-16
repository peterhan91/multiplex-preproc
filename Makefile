PY ?= python3
DATASET ?= greenwald-gbm-codex
SAMPLE ?= ZH1041_T1
CHANNEL ?= 0
MAX_WIDTH ?= 2048

.PHONY: pipeline preview panels

pipeline:
	DATASET=$(DATASET) PY=$(PY) bash run_pipeline.sh

preview:
	DATASET=$(DATASET) $(PY) scripts/render_patch_mosaic.py --sample $(SAMPLE) --channel $(CHANNEL) --max-width $(MAX_WIDTH)

panels:
	DATASET=$(DATASET) $(PY) scripts/qc_protein_panels.py --sample $(SAMPLE) --channels all
