#!/bin/bash
# MACS2 peaks + bigWigs for one group of cells (cluster / cell type / lineage) (macs2_clusters.from_signac.20240207.sh)
# Args: <group dir containing bed/> <cells_per_group table>; group = line SLURM_ARRAY_TASK_ID of the table (after header)
set -euo pipefail
source "${PIPELINE_DIR}/scripts/common.sh"
cd "$PROJECT_DIR"
activate_env "$ENV_MACS"
export LC_ALL=C

input_path=$1
input_cell_n=$2
input_filename=$(sed '1d' $input_cell_n | sed -n "${SLURM_ARRAY_TASK_ID}p" | cut -f 1)
[[ -n $input_filename ]] || die "no group at index ${SLURM_ARRAY_TASK_ID} in $input_cell_n"

output_macs=$input_path/macs
output_bam=$input_path/bam
output_bw=$input_path/bw
output_resized_bed=$input_path/resized_bed
genome_index=${GENOME_DIR}/elegans.fa.fai
output_macs_dir=$output_macs/$input_filename
mkdir -p $output_resized_bed $output_macs_dir $output_bam $output_bw

if [[ ! -s $input_path/bed/$input_filename.bed ]]; then
    log "WARNING: no fragments for group $input_filename, skipping"; exit 0
fi

awk 'BEGIN{OFS="\t";}{if ($3 - $2 < 150) print; else if ($2 - 75 < 0) print $1, 0, $2 + 75, $4, $5 "\n" $1, $3-75, $3 + 75, $4, $5; else print $1, $2 - 75, $2 + 75, $4, $5 "\n" $1, $3-75, $3 + 75, $4, $5}' $input_path/bed/$input_filename.bed | sort -k 1,1 -k2,2n > $output_resized_bed/$input_filename.resized.bed

macs2 callpeak --llocal 2000 --call-summits --bdg --SPMR --extsize 150 --shift 0 --gsize ce --keep-dup all --nomodel -n $input_filename --outdir $output_macs_dir -t $output_resized_bed/$input_filename.resized.bed
sort -k 1,1 -k2,2n $output_macs_dir/${input_filename}_treat_pileup.bdg > $output_macs_dir/${input_filename}_treat_pileup.sorted.bdg
bedGraphToBigWig $output_macs_dir/${input_filename}_treat_pileup.sorted.bdg $CHR_SIZES $output_macs_dir/${input_filename}_treat_pileup.sorted.bw
rm $output_macs_dir/${input_filename}_treat_pileup.bdg

# 100bp around summits, blacklist removed (resized again later)
slopBed -i $output_macs_dir/${input_filename}_summits.bed -b 50 -g $CHR_SIZES | intersectBed -a stdin -b $BLACKLIST -v > $output_macs_dir/${input_filename}_peaks.summits_centered.bed

# tracks: per-cell normalised and CPM
cell_n=$(awk -v var="$input_filename" '{if ($1 == var) print 1/$2}' $input_cell_n)
bedToBam -i $output_resized_bed/$input_filename.resized.bed -g $genome_index > $output_bam/$input_filename.bam
samtools index $output_bam/$input_filename.bam
bamCoverage -b $output_bam/$input_filename.bam -o $output_bw/$input_filename.cell_n_norm.bw --scaleFactor $cell_n -bs 1 -bl $BLACKLIST -p "${SLURM_CPUS_PER_TASK:-8}" --effectiveGenomeSize 100286401 --normalizeUsing None
bamCoverage -b $output_bam/$input_filename.bam -o $output_bw/$input_filename.cpm.bw -bs 1 -bl $BLACKLIST -p "${SLURM_CPUS_PER_TASK:-8}" --effectiveGenomeSize 100286401 --normalizeUsing CPM
log "group $input_filename done"
