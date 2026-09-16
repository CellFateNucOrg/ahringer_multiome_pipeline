#!/bin/bash
# Merge per-group MACS2 summits into a non-exonic peak set (pipeline blocks after "call MACS2 peaks per cluster"
# and "call MACS2 peaks per cell lineage an cell type"). Args: round1 | round2
set -euo pipefail
source "${PIPELINE_DIR}/scripts/common.sh"
cd "$PROJECT_DIR"
activate_env "$ENV_TOOLS"
export LC_ALL=C

# filter_group_peaks <dir> <narrowPeak column> <-log10 threshold> <half width> <output bed>
filter_group_peaks() {
    local dir=$1 col=$2 thr=$3 half=$4 out=$5
    cat $dir/macs/*/*.narrowPeak | awk -v c=$col -v t=$thr '{if ($c > t) print $4}' | sort > $dir/macs/all_peaks_merge.00001.names
    rm -rf $dir/intergenic_peaks; mkdir -p $dir/intergenic_peaks
    for i in $dir/macs/*/*_peaks.summits_centered.bed; do
        file_name=$(basename $i .bed)
        o=$dir/intergenic_peaks/$file_name.intergenic.bed
        # not overlapping coding/pseudogene transcripts
        intersectBed -a $i -b ${CODING_PSEUDO_BED} -v > $o
        # overlapping transcripts, but <50% overlap with exons
        intersectBed -a $i -b ${CODING_PSEUDO_BED} -u | intersectBed -a stdin -b ${CODING_PSEUDO_BED} -f 0.5 -v -split >> $o
        # remaining: keep if >=50% within a 5'UTR (or first 100bp of non-coding transcripts)
        intersectBed -a $i -b $o -v | intersectBed -a stdin -b ${UTR5_BED} -f 0.5 -u >> $o
    done
    cat $dir/intergenic_peaks/*.intergenic.bed | sort -k 4,4 | join -1 4 -2 1 - $dir/macs/all_peaks_merge.00001.names \
        | awk 'BEGIN{OFS="\t";}{print $2, $3, $4, $1}' | awk 'BEGIN{FS="\t|_"; OFS="\t";}{print $1, $2, $3, $4}' \
        | sort -k 1,1 -k2,2n | mergeBed -i stdin -c 4 -o distinct \
        | awk -v h=$half 'BEGIN{OFS="\t";}{print $1, $2 + int(($3-$2)/2) - h, $2 + int(($3-$2)/2) + h, $4}' > $out
}

case $1 in
round1)
    filter_group_peaks $ROUND1_DIR 9 6 100 $ROUND1_DIR/all_peaks_merge.intergenic.00001.tmp.bed
    intersectBed -a $ROUND1_DIR/all_peaks_merge.intergenic.00001.tmp.bed -b $BLACKLIST -v > $ROUND1_DIR/all_peaks_merge.intergenic.00001.bed
    rm $ROUND1_DIR/all_peaks_merge.intergenic.00001.tmp.bed
    log "round 1 peaks: $(wc -l < $ROUND1_DIR/all_peaks_merge.intergenic.00001.bed)"
    ;;
round2)
    for cell_annotation in cell_lineage cell_type; do
        filter_group_peaks $ROUND2_DIR/$cell_annotation 8 3 50 $ROUND2_DIR/$cell_annotation/all_peaks_merge.intergenic.00001.bed
    done
    cat $ROUND2_DIR/*/all_peaks_merge.intergenic.00001.bed | sort -k 1,1 -k2,2n | mergeBed -i stdin -c 4 -o distinct \
        | awk 'BEGIN{OFS="\t";}{print $1, $2 + int(($3-$2)/2) - 100, $2 + int(($3-$2)/2) + 100, $4}' \
        | intersectBed -a stdin -b $BLACKLIST -v > $ROUND2_DIR/all_peaks_merge.intergenic.00001.pre_IDR.bed
    awk 'BEGIN{OFS="\t";}{print $1, $2, $3, "embryo_AS_"NR}' $ROUND2_DIR/all_peaks_merge.intergenic.00001.pre_IDR.bed > $ROUND2_DIR/all_peaks.pre_IDR.bed
    log "round 2 pre-IDR peaks: $(wc -l < $ROUND2_DIR/all_peaks.pre_IDR.bed)"
    ;;
*) die "argument must be round1 or round2" ;;
esac
