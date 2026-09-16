#!/bin/bash
# WS285 transcriptome 3' extension driven by snRNA-seq coverage ("WS285 transcriptome extension" block)
set -euo pipefail
source "${PIPELINE_DIR}/scripts/common.sh"
cd "$PROJECT_DIR"
activate_env "$ENV_TOOLS"
export LC_ALL=C

base=c_elegans.PRJNA13758.${WS}.canonical_geneset
ext=$GA/extended_no_ovlp
mkdir -p $ext
summits=$BULK_PEAKS_DIR/macs/resized_fragments.all_samples_peaks.summits_centered.bed
utr_samples=$(samples_where utr_ext yes)

log "extend transcripts by ${UTR_EXT_DISTANCES} (+${UTR_EXT_BIN}bp bins)"
python scripts/gtf_3_utr_extension.20231108.py $base.filtered.gtf $UTR_EXT_BIN $UTR_EXT_DISTANCES $CHR_SIZES \
    $base.filtered.gene_transcript_id.sorted.txt $GA extended_no_ovlp

# strict extension: never extend across a scATAC accessible site
IFS=, read -r -a dists <<< "$UTR_EXT_DISTANCES"
for d in "${dists[@]}"; do
    extension=$(( d + UTR_EXT_BIN ))
    bin_size=$d
    f=$ext/$base.filtered.3prime_${bin_size}_${extension}_extended
    awk -v var="$extension" 'BEGIN{OFS="\t";}{if ($6 == "+") print $1, $3-var, $3, $4; else print $1, $2, $2 + var, $4}' $f.bed \
        | intersectBed -a stdin -b $summits -u | awk '{print $4}' | sort > $f.to_remove
    sort -k 4,4 $f.bed | join -v 1 -1 4 -2 1 - $f.to_remove | awk 'BEGIN{OFS="\t";}{print $2, $3, $4, $1, $5, $6}' | grep -v MtDNA > $f.no_atac_peaks.bed
done

cat $ext/$base.filtered.3prime_*_extended.no_atac_peaks.bed | sortBed -i stdin -faidx $CHR_SIZES > $ext/$base.filtered.3prime_extended.all.whole_transcript.bed

# coverage of deduplicated snRNA-seq reads: extended bins (blacklist removed) and annotated transcripts (whole locus, bed6)
grep -v MtDNA $GA/$base.filtered.bed | cut -f 1,2,3,4,5,6 | sort -k 1,1 -k2,2n > $GA/$base.filtered.no_MtDNA.bed6
ext_cov=(); orig_cov=()
for id in $utr_samples; do
    bam=data/sc/rnaseq/$id/star_${WS}/star_${WS}_Aligned.sortedByCoord.out.split.merge.dedup.bam
    log "coverage $id"
    intersectBed -a $ext/$base.filtered.3prime_extended.all.whole_transcript.bed -b $BLACKLIST -v | sort -k 1,1 -k2,2n \
        | coverageBed -a stdin -b $bam -counts -s -g $CHR_SIZES -sorted > $ext/$base.filtered.3prime_extended.all.whole_transcript.$id.coverage
    coverageBed -a $GA/$base.filtered.no_MtDNA.bed6 -b $bam -counts -s -g $CHR_SIZES -sorted > $ext/$base.filtered.$id.coverage
    ext_cov+=("$ext/$base.filtered.3prime_extended.all.whole_transcript.$id.coverage")
    orig_cov+=("$ext/$base.filtered.$id.coverage")
done

sum_cov() { paste "$@" | awk -v n=$# 'BEGIN{OFS="\t";}{s=0; for (i=0; i<n; i++) s+=$(7+7*i); print $1, $2, $3, $4, $5, $6, s}'; }
sum_cov "${ext_cov[@]}" > $ext/$base.filtered.3prime_extended.all.whole_transcript.all_samples.coverage
sum_cov "${orig_cov[@]}" > $ext/$base.filtered.all_samples.coverage

log "select extensions"
python scripts/gtf_3_utr_extension.bin_selection.20231108.py $ext/$base.filtered.all_samples.coverage \
    $ext/$base.filtered.3prime_extended.all.whole_transcript.all_samples.coverage $UTR_EXT_BIN $UTR_EXT_DISTANCES \
    $ext/$base.filtered.3prime_extensions.txt

python scripts/gtf_3_utr_extension.rnaseq_extended.20231108.py $base.filtered.bed $base.filtered.3prime_extensions.txt \
    $UTR_EXT_BIN $UTR_EXT_DISTANCES $base.filtered.gene_transcript_id.sorted.txt $GA extended_no_ovlp
python scripts/gtf_3_utr_extension.rnaseq_extended.20231108.py $base.bed $base.filtered.3prime_extensions.txt \
    $UTR_EXT_BIN $UTR_EXT_DISTANCES $base.gene_transcript_id.sorted.txt $GA extended_no_ovlp

# coding + pseudogene transcripts of the extended annotation
sort -k 1,1 ${CANON}.gene_transcript_id.sorted.txt | join - ${CANON}.gene_type.txt | awk '{if ($3 == "protein_coding" || $3 == "pseudogene") print $2}' | sort \
    | join -1 1 -2 4 - <(sort -k 4,4 ${CANON}.extended_no_ovlp.3prime_extended.rnaseq_corrected.bed) \
    | awk 'BEGIN{OFS="\t";}{print $2, $3, $4, $1, $5, $6, $7, $8, $9, $10, $11, $12}' > ${CODING_PSEUDO_BED}

[[ -s ${CANON}.extended_no_ovlp.3prime_extended.rnaseq_corrected.gtf ]] || die "extended GTF was not produced"
log "extended GTF: ${CANON}.extended_no_ovlp.3prime_extended.rnaseq_corrected.gtf ($(awk '$2 != 0' $ext/$base.filtered.3prime_extensions.txt | wc -l) transcripts extended)"
