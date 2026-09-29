#!/bin/bash
# Split one sample's STARsolo BAM (extended annotation) by cell type (bam_barcode_split.nosuffix.sh)
set -euo pipefail
source "${PIPELINE_DIR}/scripts/common.sh"
cd "$PROJECT_DIR"
activate_env "$ENV_TOOLS"

id=$(sample_by_index "${SLURM_ARRAY_TASK_ID}")
bam=data/sc/rnaseq/$id/star_${ANNOTATION}/star_${ANNOTATION}_Aligned.sortedByCoord.out.bam
barcodes=tracks/rna_barcodes/$id.barcodes.txt
[[ -s $barcodes ]] || { log "no annotated cells for $id"; exit 0; }
out=tracks/rna_bam_split/$id
mkdir -p $out

# skip if all expected cell-type BAMs already exist
n_expected=$(awk '{print $2}' $barcodes | sort -u | wc -l)
n_actual=$(ls $out/*.bam 2>/dev/null | wc -l)
if [[ $n_actual -ge $n_expected && $n_expected -gt 0 ]]; then
    log "$id: already split into $n_actual cell types, skipping"
    exit 0
fi

[[ -s $bam.bai ]] || samtools index -@ "${SLURM_CPUS_PER_TASK:-12}" $bam
sinto filterbarcodes -b $bam -c $barcodes --outdir $out -p "${SLURM_CPUS_PER_TASK:-12}"
log "$id split into $(ls $out/*.bam | wc -l) cell types"
