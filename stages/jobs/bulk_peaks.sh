#!/bin/bash
# MACS2 peaks on pooled scATAC fragments of the utr_ext samples (macs2_cellranger_fragments_all.20231103.sh)
set -euo pipefail
source "${PIPELINE_DIR}/scripts/common.sh"
cd "$PROJECT_DIR"
activate_env "$ENV_MACS"
export LC_ALL=C

output_dir=$BULK_PEAKS_DIR
mkdir -p $output_dir/macs $output_dir/tmp
input_files=()
for id in $(samples_where utr_ext yes); do input_files+=("data/sc/cellranger_arc/$id/outs/atac_fragments.tsv.gz"); done
log "fragments: ${input_files[*]}"

# split each fragment >=150bp into its two Tn5 cut sites, each extended to 150bp
zcat "${input_files[@]}" | grep -v "#" | awk 'BEGIN{OFS="\t";}{if ($3 - $2 < 150) print; else if ($2 - 75 < 0) print $1, 0, $2 + 75, $4, $5 "\n" $1, $3-75, $3 + 75, $4, $5; else print $1, $2 - 75, $2 + 75, $4, $5 "\n" $1, $3-75, $3 + 75, $4, $5}' \
    | sort -S 60% -T $output_dir/tmp -k 1,1 -k2,2n > $output_dir/resized_fragments.all.bed

macs2 callpeak --call-summits --bdg --SPMR --extsize 150 --shift 0 --gsize ce --keep-dup all --nomodel -n resized_fragments.all_samples --outdir $output_dir/macs -t $output_dir/resized_fragments.all.bed
sort -k 1,1 -k2,2n $output_dir/macs/resized_fragments.all_samples_treat_pileup.bdg > $output_dir/macs/resized_fragments.all_samples_treat_pileup.sorted.bdg
awk 'BEGIN{OFS="\t";}{print $1, 1, $2}' $CHR_SIZES | intersectBed -a $output_dir/macs/resized_fragments.all_samples_treat_pileup.sorted.bdg -b stdin -u -f 1 > $output_dir/macs/resized_fragments.all_samples_treat_pileup.sorted.filtered.bdg
bedGraphToBigWig $output_dir/macs/resized_fragments.all_samples_treat_pileup.sorted.filtered.bdg $CHR_SIZES $output_dir/macs/resized_fragments.all_samples_treat_pileup.sorted.bw
rm $output_dir/macs/resized_fragments.all_samples_treat_pileup.bdg $output_dir/macs/resized_fragments.all_samples_treat_pileup.sorted.bdg

# 200bp peaks centred on summits
slopBed -i $output_dir/macs/resized_fragments.all_samples_summits.bed -b 100 -g $CHR_SIZES > $output_dir/macs/resized_fragments.all_samples_peaks.summits_centered.bed
cut -f 1,2,3,4,5,6 $output_dir/macs/resized_fragments.all_samples_peaks.narrowPeak > $output_dir/macs/resized_fragments.all_samples_peaks.bed
rm -rf $output_dir/tmp
log "bulk peaks done: $(wc -l < $output_dir/macs/resized_fragments.all_samples_peaks.bed) peaks"
