#!/bin/bash
# Fragment file preparation for one sample:
#  1. remove the "-1" suffix from barcodes (atac_fragments_tsv_parser.20230502.py)
#  2. prefix barcodes with the sample id (atac_fragment_sample_barcode.20231107.sh)
set -euo pipefail
source "${PIPELINE_DIR}/scripts/common.sh"
cd "$PROJECT_DIR"
activate_env "$ENV_TOOLS"

sample_id=$(sample_by_index "${SLURM_ARRAY_TASK_ID}")
outs=data/sc/cellranger_arc/$sample_id/outs

no_dash=$outs/atac_fragments_no_dash.tsv.gz
corrected=$outs/atac_fragments_barcode_corrected.tsv.gz

if [[ -s $corrected && -s ${corrected}.tbi ]]; then
    log "$sample_id: fragment files exist, skipping"; exit 0
fi

zcat $outs/atac_fragments.tsv.gz | awk 'BEGIN{OFS="\t"} /^#/{print; next} {sub(/-.*$/, "", $4); print}' | bgzip -f > $no_dash
tabix -f -p bed $no_dash

zcat $no_dash | awk -v var="$sample_id" 'BEGIN{OFS="\t";} /^#/{print; next} {print $1, $2, $3, var "_" $4, $5}' | bgzip -f > $corrected
tabix -f -p bed $corrected
log "$sample_id: fragment files written"
