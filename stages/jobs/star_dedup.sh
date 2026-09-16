#!/bin/bash
# UMI deduplication of the STARsolo BAM, split into 20 barcode groups to limit memory
# (bam_barcode_split.20231031.sh, STAR_dedup.20231031.sh, merge_dedup_bam.20231031.sh). Args: <index name>
# Only for samples with utr_ext=yes (used for 3' UTR extension coverage).
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

split_dir=data/sc/rnaseq/$id/star_$idx/star_${idx}_Aligned.barcode_split
mkdir -p $split_dir
samtools index -@ $threads $bam

# assign whitelist barcodes to 20 groups (00..19); the grouping only affects memory use, not the result
total=$(wc -l < whitelist/737K-arc-v1.txt)
awk -v n=$total 'BEGIN{OFS="\t"; chunk=int((n+19)/20)} {printf "%s\t%02d\n", $1, int((NR-1)/chunk)}' whitelist/737K-arc-v1.txt > $split_dir/barcode_groups.txt
sinto filterbarcodes -b $bam -c $split_dir/barcode_groups.txt --outdir $split_dir -p $threads

ls $split_dir/[0-9][0-9].bam | xargs -P "${DEDUP_PARALLEL}" -I{} bash -c '
    f="$1"; samtools index "$f"
    umi_tools dedup --stdin="$f" --log="${f%.bam}.dedup.log" --error="${f%.bam}.dedup.err" \
        --extract-umi-method=tag --umi-tag=UR --cell-tag=CR --per-cell > "${f%.bam}.dedup.bam"' _ {}

samtools merge -f -h $bam -@ $threads $final $split_dir/[0-9][0-9].dedup.bam
samtools index -@ $threads $final
log "$id: deduplicated BAM $final"
