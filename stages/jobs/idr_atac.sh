#!/bin/bash
# IDR between pseudo-replicates for one cell type (idr_atac.sh + idr_vs_all_macs_peaks_atac_comparison.sh)
# Args: <cells_per_cell_type table>; cell type = line SLURM_ARRAY_TASK_ID (after header)
set -euo pipefail
source "${PIPELINE_DIR}/scripts/common.sh"
cd "$PROJECT_DIR"
activate_env "$ENV_MACS"

input_path=$ROUND2_DIR/cell_type/macs_peak_frag_counts
input_filename=$(sed '1d' $1 | sed -n "${SLURM_ARRAY_TASK_ID}p" | cut -f 1)
if ! grep -qx -- "$input_filename" $input_path/cell_types.txt; then
    log "$input_filename not present in both pseudo-replicates, skipping"; exit 0
fi
mkdir -p $input_path/IDR $input_path/IDR_logs

if ! idr --samples $input_path/${input_filename}_rep_1.bed $input_path/${input_filename}_rep_2.bed --input-file-type narrowPeak \
        --output-file $input_path/IDR/${input_filename}_idr.txt --output-file-type bed \
        --log-output-file $input_path/IDR_logs/${input_filename}_idr.log --plot --rank signal.value; then
    log "WARNING: IDR failed for $input_filename (see IDR_logs); no peaks retained for this cell type"
    : > $input_path/IDR/${input_filename}_idr.txt
fi

awk -v var="$input_filename" 'BEGIN{OFS="\t";}{if($5 >= 830) print $1, $2, $3, var}' $input_path/IDR/${input_filename}_idr.txt > $input_path/IDR/${input_filename}_idr.0.01.bed
awk -v var="$input_filename" 'BEGIN{OFS="\t";}{if($5 >= 540) print $1, $2, $3, var}' $input_path/IDR/${input_filename}_idr.txt > $input_path/IDR/${input_filename}_idr.0.05.bed

bw=$ROUND2_DIR/cell_type/bw/$input_filename.cpm.bw
if [[ ${IDR_QC_HEATMAPS} == "yes" && -s $bw && -s $input_path/IDR/${input_filename}_idr.${IDR_FINAL_SET}.bed ]]; then
    out=$input_path/IDR_vs_MACS; mkdir -p $out
    q=$out/$input_filename.idr.${IDR_FINAL_SET}.macs
    intersectBed -a $input_path/IDR/${input_filename}_idr.${IDR_FINAL_SET}.bed -b $ROUND2_DIR/all_peaks.pre_IDR.bed -u | cut -f 1,2,3 > $q.bed
    echo "#IDR" >> $q.bed
    intersectBed -b $input_path/IDR/${input_filename}_idr.${IDR_FINAL_SET}.bed -a $ROUND2_DIR/all_peaks.pre_IDR.bed -v | cut -f 1,2,3 >> $q.bed
    echo "#no_IDR" >> $q.bed
    computeMatrix scale-regions -R $q.bed -S $bw -o $q.mtx.gz -a 500 -b 500 -m 200 --binSize 10 -p "${SLURM_CPUS_PER_TASK:-8}" -bl $BLACKLIST
    plotHeatmap -m $q.mtx.gz -o $q.pdf --heatmapHeight 25 --heatmapWidth 5 --colorMap jet --whatToShow 'heatmap and colorbar' --zMin 0 --zMax 20
fi
log "$input_filename: $(wc -l < $input_path/IDR/${input_filename}_idr.${IDR_FINAL_SET}.bed) peaks at IDR ${IDR_FINAL_SET}"
