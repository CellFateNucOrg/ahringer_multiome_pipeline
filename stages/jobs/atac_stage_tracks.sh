#!/bin/bash
# Average ATAC CPM tracks of cell types belonging to the same embryo stage (combine_bw_tracks_stage.20240916.sh),
# using the 'stage' metadata written by seurat_integrate_cluster_peaks.R instead of hard-coded cell type lists
set -euo pipefail
source "${PIPELINE_DIR}/scripts/common.sh"
cd "$PROJECT_DIR"
activate_env "$ENV_MACS"

out=$ROUND2_DIR/embryo_stage_tracks
mkdir -p $out
for stage in $(cut -f 1 tracks/stage_cell_types.tsv | sort -u); do
    bws=()
    for ct in $(awk -F'\t' -v s="$stage" '$1 == s {print $2}' tracks/stage_cell_types.tsv); do
        f=$ROUND2_DIR/cell_type/bw/$ct.cpm.bw
        [[ -s $f ]] && bws+=("$f")
    done
    if [[ ${#bws[@]} -eq 0 ]]; then
        log "stage $stage: no tracks"
    elif [[ ${#bws[@]} -eq 1 ]]; then
        cp "${bws[0]}" $out/$stage.cpm.bw
    else
        bigwigAverage --bigwigs "${bws[@]}" -bs 1 -p "${SLURM_CPUS_PER_TASK:-12}" -o $out/$stage.cpm.bw
    fi
done
log "stage tracks in $out"
