#!/bin/bash
# Merge per-sample BAMs of one cell type, deduplicate UMIs and write strand-specific CPM bigWigs
# (merge_replicate_bam.20240502.sh + dedup_launcher.clusters.20240226.sh). Args: <cells_per_cell_type table>
set -euo pipefail
source "${PIPELINE_DIR}/scripts/common.sh"
cd "$PROJECT_DIR"
activate_env "$ENV_TOOLS"

sample_name=$(sed '1d' $1 | sed -n "${SLURM_ARRAY_TASK_ID}p" | cut -f 1)
alignment_path=tracks/rna
threads=${SLURM_CPUS_PER_TASK:-8}
mkdir -p $alignment_path/bam $alignment_path/dedup_err $alignment_path/dedup_log $alignment_path/bw/unique $alignment_path/dedup_bam

# skip if all final bigWigs already exist
d=$alignment_path/dedup_bam/$sample_name.dedup.bam
if [[ -s $d && \
      -s $alignment_path/bw/$sample_name.dedup.rev.bw && \
      -s $alignment_path/bw/$sample_name.dedup.for.bw && \
      -s $alignment_path/bw/unique/$sample_name.unique.dedup.rev.bw && \
      -s $alignment_path/bw/unique/$sample_name.unique.dedup.for.bw ]]; then
    log "$sample_name: tracks already exist, skipping"
    exit 0
fi

mapfile -t inputs < <(ls tracks/rna_bam_split/*/"$sample_name".bam 2>/dev/null || true)
[[ ${#inputs[@]} -gt 0 ]] || { log "no BAM for $sample_name, skipping"; exit 0; }

# merge if needed
merged=$alignment_path/bam/$sample_name.bam
if [[ ! -s $merged ]]; then
    samtools merge -f -@ $threads $merged "${inputs[@]}"
    samtools index -@ $threads $merged
fi

# dedup if needed
if [[ ! -s $d ]]; then
    umi_tools dedup --stdin=$merged --log=$alignment_path/dedup_log/$sample_name.dedup.log \
        --error=$alignment_path/dedup_err/$sample_name.dedup.err \
        --extract-umi-method=tag --umi-tag=UR --cell-tag=CR --per-cell > $d
    samtools index $d
fi

# bigWigs
[[ -s $alignment_path/bw/$sample_name.dedup.rev.bw ]] || \
    bamCoverage -b $d -o $alignment_path/bw/$sample_name.dedup.rev.bw --filterRNAstrand forward -bs 1 -p $threads \
        --effectiveGenomeSize 100286401 --normalizeUsing CPM --blackListFileName=$BLACKLIST
[[ -s $alignment_path/bw/unique/$sample_name.unique.dedup.rev.bw ]] || \
    bamCoverage -b $d -o $alignment_path/bw/unique/$sample_name.unique.dedup.rev.bw --filterRNAstrand forward -bs 1 -p $threads \
        --effectiveGenomeSize 100286401 --normalizeUsing CPM --samFlagExclude 256 --blackListFileName=$BLACKLIST
[[ -s $alignment_path/bw/$sample_name.dedup.for.bw ]] || \
    bamCoverage -b $d -o $alignment_path/bw/$sample_name.dedup.for.bw --filterRNAstrand reverse -bs 1 -p $threads \
        --effectiveGenomeSize 100286401 --normalizeUsing CPM --blackListFileName=$BLACKLIST
[[ -s $alignment_path/bw/unique/$sample_name.unique.dedup.for.bw ]] || \
    bamCoverage -b $d -o $alignment_path/bw/unique/$sample_name.unique.dedup.for.bw --filterRNAstrand reverse -bs 1 -p $threads \
        --effectiveGenomeSize 100286401 --normalizeUsing CPM --samFlagExclude 256 --blackListFileName=$BLACKLIST
log "$sample_name RNA tracks done"
