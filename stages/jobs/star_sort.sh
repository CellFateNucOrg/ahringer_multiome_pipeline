#!/bin/bash
# Sort the unsorted STARsolo BAM with samtools, using local scratch for temp files.
# Args: <index name>
set -euo pipefail
source "${PIPELINE_DIR}/scripts/common.sh"
cd "$PROJECT_DIR"
activate_env "$ENV_TOOLS"

idx=$1
id=$(sample_by_index "${SLURM_ARRAY_TASK_ID}")
outdir=data/sc/rnaseq/$id/star_$idx
prefix=$outdir/star_${idx}_

unsorted=${prefix}Aligned.out.bam
sorted=${prefix}Aligned.sortedByCoord.out.bam

if [[ -s $sorted ]]; then
    log "$id: sorted BAM exists, skipping sort"
    exit 0
fi

[[ -s $unsorted ]] || die "$id: unsorted BAM not found: $unsorted (did starsolo.sh complete?)"

threads=${SLURM_CPUS_PER_TASK:-8}
# allocate 75% of job RAM across threads for sorting (in bytes per thread)
mem_per_thread=$(( ${SLURM_MEM_PER_NODE:-131072} * 1024 * 1024 * 3 / 4 / threads ))

if [[ -n "${SCRATCH_DIR:-}" ]]; then
    scratch="${SCRATCH_DIR}/star_sort_${SLURM_JOB_ID:-0}_${SLURM_ARRAY_TASK_ID:-0}"
    mkdir -p "$scratch"
    trap 'rm -rf "$scratch"' EXIT
    sort_tmp_opt="-T $scratch/sort"
    log "$id: sorting BAM (threads=$threads, mem_per_thread=$(( mem_per_thread / 1024 / 1024 ))MB, scratch=$scratch)"
else
    sort_tmp_opt=""
    log "$id: sorting BAM (threads=$threads, mem_per_thread=$(( mem_per_thread / 1024 / 1024 ))MB, no scratch dir)"
fi

# shellcheck disable=SC2086
samtools sort -@ "$threads" -m "${mem_per_thread}" $sort_tmp_opt \
    -o "$sorted" "$unsorted"
samtools index -@ "$threads" "$sorted"

rm -f "$unsorted"
log "$id: BAM sort done → $sorted"

# Generate whole-sample bigWigs (forward and reverse strand, CPM-normalised)
bw_dir=$outdir/bw
mkdir -p "$bw_dir"
fwd_bw=$bw_dir/${id}.star_${idx}_fwd.bw
rev_bw=$bw_dir/${id}.star_${idx}_rev.bw

if [[ -s $fwd_bw && -s $rev_bw ]]; then
    log "$id: bigWigs exist, skipping"
else
    log "$id: generating bigWigs"
    # Note: deeptools --filterRNAstrand forward captures reverse-strand reads (library convention)
    bamCoverage -b "$sorted" -o "$fwd_bw" --filterRNAstrand reverse \
        -bs 1 -p "$threads" --effectiveGenomeSize 100286401 \
        --normalizeUsing CPM --blackListFileName="$BLACKLIST"
    bamCoverage -b "$sorted" -o "$rev_bw" --filterRNAstrand forward \
        -bs 1 -p "$threads" --effectiveGenomeSize 100286401 \
        --normalizeUsing CPM --blackListFileName="$BLACKLIST"
    log "$id: bigWigs done → $bw_dir"
fi
