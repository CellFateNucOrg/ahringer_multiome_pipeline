#!/bin/bash
# UMI deduplication of the STARsolo BAM, split into 20 barcode groups to limit memory
# (bam_barcode_split.20231031.sh, STAR_dedup.20231031.sh, merge_dedup_bam.20231031.sh). Args: <index name>
# Only for samples with utr_ext=yes (used for 3' UTR extension coverage).
#
# Scratch strategy: split BAMs (large, temporary) go to node-local scratch for fast I/O.
# Dedup outputs go directly to NFS so they persist as resume checkpoints if the job is
# killed — on the next run sinto re-splits, but already-deduped groups are skipped.
set -euo pipefail
source "${PIPELINE_DIR}/scripts/common.sh"
cd "$PROJECT_DIR"
activate_env "$ENV_TOOLS"

idx=$1
id=$(sample_by_index "${SLURM_ARRAY_TASK_ID}")
[[ $(sample_get "$id" utr_ext) == "yes" ]] || { log "$id: utr_ext != yes, nothing to do"; exit 0; }

bam=data/sc/rnaseq/$id/star_$idx/star_${idx}_Aligned.sortedByCoord.out.bam
final=data/sc/rnaseq/$id/star_$idx/star_${idx}_Aligned.sortedByCoord.out.split.merge.dedup.bam
[[ -s $final ]] && { log "$id: $final exists, skipping"; exit 0; }
threads=${SLURM_CPUS_PER_TASK:-8}

# NFS dir: holds barcode_groups.txt, dedup BAMs, logs — persists across runs
split_dir=data/sc/rnaseq/$id/star_$idx/star_${idx}_Aligned.barcode_split
mkdir -p $split_dir

# ── Scratch setup ────────────────────────────────────────────────────────────
if [[ -n "${SCRATCH_DIR:-}" ]]; then
    scratch="${SCRATCH_DIR}/star_dedup_${id}_${idx}"
    if [[ -d "$scratch" ]]; then
        log "$id: removing stale scratch dir $scratch from a previous (killed?) run"
        rm -rf "$scratch"
    fi
    mkdir -p "$scratch"
    trap 'log "$id: dedup failed or interrupted; cleaning scratch $scratch"; rm -rf "$scratch"' ERR EXIT
    split_scratch="$scratch/barcode_split"
    log "$id: split BAMs → scratch $split_scratch; dedup BAMs → NFS $split_dir"
else
    split_scratch="$split_dir"
    log "$id: no scratch dir; split BAMs and dedup BAMs both on NFS"
fi
mkdir -p "$split_scratch"

samtools index -@ $threads $bam

# assign whitelist barcodes to 20 groups (00..19); the grouping only affects memory use, not the result
total=$(wc -l < whitelist/737K-arc-v1.txt)
awk -v n=$total 'BEGIN{OFS="\t"; chunk=int((n+19)/20)} {printf "%s\t%02d\n", $1, int((NR-1)/chunk)}' \
    whitelist/737K-arc-v1.txt > $split_dir/barcode_groups.txt

# ── sinto: split BAM into 20 barcode groups → scratch ───────────────────────
# Skip sinto if all 20 scratch split BAMs already exist (within-job resume after ERR).
# On a fresh run or after a SIGKILL (scratch wiped above), sinto always re-runs;
# the per-group dedup skip below avoids re-doing completed groups.
n_splits=$(ls "$split_scratch"/[0-9][0-9].bam 2>/dev/null | wc -l)
if [[ $n_splits -eq 20 ]]; then
    log "$id: all 20 barcode-split BAMs present in scratch, skipping sinto"
else
    log "$id: running sinto filterbarcodes"
    sinto filterbarcodes -b $bam -c $split_dir/barcode_groups.txt --barcodetag CR \
        --outdir "$split_scratch" -p $threads
fi

# ── per-BAM dedup ─────────────────────────────────────────────────────────────
# Resume checkpoint is the NFS dedup BAM. Each split BAM is removed from scratch
# immediately after its dedup BAM lands on NFS, freeing space progressively.
export split_dir split_scratch
ls "$split_scratch"/[0-9][0-9].bam | xargs -P "${DEDUP_PARALLEL}" -I{} bash -c '
    scratch_bam="$1"
    group=$(basename "${scratch_bam%.bam}")
    nfs_dedup="$split_dir/$group.dedup.bam"
    if [[ -s "$nfs_dedup" ]]; then
        echo "$(date "+%Y-%m-%d %H:%M:%S") $group.dedup.bam already on NFS, skipping"
        rm -f "$scratch_bam" "${scratch_bam}.bai"
    else
        samtools index "$scratch_bam"
        umi_tools dedup --stdin="$scratch_bam" \
            --log="$split_dir/$group.dedup.log" --error="$split_dir/$group.dedup.err" \
            --extract-umi-method=tag --umi-tag=UR --cell-tag=CR --per-cell > "$nfs_dedup"
        rm -f "$scratch_bam" "${scratch_bam}.bai"
        echo "$(date "+%Y-%m-%d %H:%M:%S") $group.dedup.bam written to NFS; scratch split BAM removed"
    fi' _ {}

samtools merge -f -h $bam -@ $threads $final $split_dir/[0-9][0-9].dedup.bam
samtools index -@ $threads $final
trap - ERR EXIT
[[ -n "${SCRATCH_DIR:-}" ]] && rm -rf "$scratch"
log "$id: deduplicated BAM $final"
