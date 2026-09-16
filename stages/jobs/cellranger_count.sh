#!/bin/bash
# cellranger-arc count for one sample (array task index = sample row in samples.tsv)
set -euo pipefail
source "${PIPELINE_DIR}/scripts/common.sh"
cd "$PROJECT_DIR"

id=$(sample_by_index "${SLURM_ARRAY_TASK_ID}")
outs=data/sc/cellranger_arc/$id/outs
if [[ -s $outs/atac_fragments.tsv.gz && -s $outs/raw_feature_bc_matrix.h5 ]]; then
    log "$id: cellranger-arc output exists, skipping"; exit 0
fi

lib=cellranger_arc_library_files/$id.cr_arc_library_file.csv
mkdir -p cellranger_arc_library_files
echo "fastqs,sample,library_type" > $lib
IFS=, read -r -a gex <<< "$(sample_get "$id" gex_fastq_prefix)"
for p in "${gex[@]}"; do echo "$(sample_get "$id" gex_fastq_dir),$p,Gene Expression" >> $lib; done
IFS=, read -r -a atac <<< "$(sample_get "$id" atac_fastq_prefix)"
for p in "${atac[@]}"; do echo "$(sample_get "$id" atac_fastq_dir),$p,Chromatin Accessibility" >> $lib; done
cat $lib

mem_gb=$(( ${SLURM_MEM_PER_NODE:-64000} / 1024 - 4 ))
cd data/sc/cellranger_arc
"${CELLRANGER_ARC_DIR}/cellranger-arc" count --id=$id --reference="$PROJECT_DIR/$CR_REF" \
    --libraries="$PROJECT_DIR/$lib" --localcores="${SLURM_CPUS_PER_TASK:-12}" --localmem=$mem_gb --create-bam true
log "$id: cellranger-arc count done"
