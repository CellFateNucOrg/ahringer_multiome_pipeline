#!/bin/bash
# Final accessible-site set = union of IDR-reproducible peaks of all cell types ("define final peakset")
set -euo pipefail
source "${PIPELINE_DIR}/scripts/common.sh"
cd "$PROJECT_DIR"
activate_env "$ENV_TOOLS"
export LC_ALL=C

idr_dir=$ROUND2_DIR/cell_type/macs_peak_frag_counts/IDR
cat $idr_dir/*_idr.${IDR_FINAL_SET}.bed | sort -k 1,1 -k2,2n | mergeBed -i stdin -c 4 -o distinct \
    | awk 'BEGIN{OFS="\t";}{print $1, $2 + int(($3-$2)/2) - 100, $2 + int(($3-$2)/2) + 100, $4}' \
    | intersectBed -a stdin -b $BLACKLIST -v > $ROUND2_DIR/all_peaks_merge.intergenic.00001.bed
awk 'BEGIN{OFS="\t";}{print $1, $2, $3, "embryo_AS_"NR}' $ROUND2_DIR/all_peaks_merge.intergenic.00001.bed > $ROUND2_DIR/all_peaks.bed
fastaFromBed -bed $ROUND2_DIR/all_peaks.bed -fi ${GENOME_DIR}/elegans.fa -fo $ROUND2_DIR/all_peaks.fa

# accessible sites per cell type
mkdir -p $ROUND2_DIR/cell_type/final_peakset
for i in $idr_dir/*_idr.${IDR_FINAL_SET}.bed; do
    file_name=$(basename $i _idr.${IDR_FINAL_SET}.bed)
    intersectBed -a $ROUND2_DIR/all_peaks.bed -b $i -u > $ROUND2_DIR/cell_type/final_peakset/$file_name.peaks.bed
done
log "final peak set: $(wc -l < $ROUND2_DIR/all_peaks.bed) accessible sites"
