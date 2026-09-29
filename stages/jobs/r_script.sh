#!/bin/bash
# Generic R job: r_script.sh <script.R> [args...]
# The literal argument __SAMPLE__ is replaced by the sample id of the Slurm array task.
set -euo pipefail
source "${PIPELINE_DIR}/scripts/common.sh"
cd "$PROJECT_DIR"
activate_env "$ENV_R"

export SAMPLES_TSV MT_PCT_MAX MAX_RNA_UMI SOUPX_QUANTILE BATCH_FOLD LATE_STAGE_GROUPS LATE_MAX_FRACTION \
       MIN_WNN_FEATURES UNASSIGN_CELL_TYPES WS ANNOTATION MARKER_LIST_FILE MARKER_TAXONOMY_FILE

script=$1; shift
args=()
for a in "$@"; do
    [[ $a == "__SAMPLE__" ]] && a=$(sample_by_index "${SLURM_ARRAY_TASK_ID}")
    args+=("$a")
done
export R_LIBS_USER=/nonexistent
export R_LIBS=
log "Rscript scripts/$script ${args[*]:-}"
Rscript --no-environ --no-site-file --no-init-file "scripts/$script" ${args[@]+"${args[@]}"}
log "done"
