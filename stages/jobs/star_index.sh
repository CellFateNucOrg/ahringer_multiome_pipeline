#!/bin/bash
# STAR genome index (STARsolo_indexing.sh). Args: <gtf> <index name>
set -euo pipefail
source "${PIPELINE_DIR}/scripts/common.sh"
cd "$PROJECT_DIR"
activate_env "$ENV_TOOLS"

input_gtf=$1
genome_index=$2
mkdir -p star_idx/$genome_index

[[ -s star_idx/$genome_index/SA ]] && { log "STAR index star_idx/$genome_index already exists, skipping"; exit 0; }

STAR --runMode genomeGenerate --genomeSAindexNbases 12 --sjdbGTFfile $input_gtf --runThreadN "${SLURM_CPUS_PER_TASK:-10}" \
     --genomeDir star_idx/$genome_index --genomeFastaFiles ${GENOME_DIR}/elegans.fa --outFileNamePrefix star_idx/$genome_index/
log "STAR index star_idx/$genome_index done"
