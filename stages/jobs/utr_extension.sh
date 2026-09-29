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

# ── Overall skip ─────────────────────────────────────────────────────────────
if [[ -s ${CANON}.extended_no_ovlp.3prime_extended.rnaseq_corrected.gtf && -s ${CODING_PSEUDO_BED} ]]; then
    log "UTR extension already complete, skipping"
    exit 0
fi

# ── Step 1: compute extension bins ──────────────────────────────────────────
if [[ -s $ext/$base.filtered.3prime_extended.all.whole_transcript.bed ]]; then
    log "extension bins already computed, skipping step 1"
else
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
fi

# ── Step 2: coverage of deduplicated snRNA-seq reads ────────────────────────
grep -v MtDNA $GA/$base.filtered.bed | cut -f 1,2,3,4,5,6 | sort -k 1,1 -k2,2n > $GA/$base.filtered.no_MtDNA.bed6
ext_cov=(); orig_cov=()
for id in $utr_samples; do
    bam=data/sc/rnaseq/$id/star_${WS}/star_${WS}_Aligned.sortedByCoord.out.split.merge.dedup.bam
    ec=$ext/$base.filtered.3prime_extended.all.whole_transcript.$id.coverage
    oc=$ext/$base.filtered.$id.coverage
    if [[ -s $ec && -s $oc ]]; then
        log "coverage $id already done, skipping"
    else
        log "coverage $id"
        intersectBed -a $ext/$base.filtered.3prime_extended.all.whole_transcript.bed -b $BLACKLIST -v | sort -k 1,1 -k2,2n \
            | coverageBed -a stdin -b $bam -counts -s -g $CHR_SIZES -sorted > $ec
        coverageBed -a $GA/$base.filtered.no_MtDNA.bed6 -b $bam -counts -s -g $CHR_SIZES -sorted > $oc
    fi
    ext_cov+=("$ec")
    orig_cov+=("$oc")
done

sum_cov() { paste "$@" | awk -v n=$# 'BEGIN{OFS="\t";}{s=0; for (i=0; i<n; i++) s+=$(7+7*i); print $1, $2, $3, $4, $5, $6, s}'; }
[[ -s $ext/$base.filtered.3prime_extended.all.whole_transcript.all_samples.coverage ]] || \
    sum_cov "${ext_cov[@]}" > $ext/$base.filtered.3prime_extended.all.whole_transcript.all_samples.coverage
[[ -s $ext/$base.filtered.all_samples.coverage ]] || \
    sum_cov "${orig_cov[@]}" > $ext/$base.filtered.all_samples.coverage

# ── Step 3: select extensions ─────────────────────────────────────────────────
if [[ ! -s $ext/$base.filtered.3prime_extensions.txt ]]; then
    log "select extensions"
    python scripts/gtf_3_utr_extension.bin_selection.20231108.py $ext/$base.filtered.all_samples.coverage \
        $ext/$base.filtered.3prime_extended.all.whole_transcript.all_samples.coverage $UTR_EXT_BIN $UTR_EXT_DISTANCES \
        $ext/$base.filtered.3prime_extensions.txt
fi

# ── Step 4: rnaseq-corrected extended annotation ──────────────────────────────
if [[ ! -s ${CANON}.extended_no_ovlp.3prime_extended.rnaseq_corrected.gtf ]]; then
    python scripts/gtf_3_utr_extension.rnaseq_extended.20231108.py $base.filtered.bed $base.filtered.3prime_extensions.txt \
        $UTR_EXT_BIN $UTR_EXT_DISTANCES $base.filtered.gene_transcript_id.sorted.txt $GA extended_no_ovlp
    python scripts/gtf_3_utr_extension.rnaseq_extended.20231108.py $base.bed $base.filtered.3prime_extensions.txt \
        $UTR_EXT_BIN $UTR_EXT_DISTANCES $base.gene_transcript_id.sorted.txt $GA extended_no_ovlp
fi

# ── Step 5: coding + pseudogene transcripts of the extended annotation ────────
if [[ ! -s ${CODING_PSEUDO_BED} ]]; then
    sort -k 1,1 ${CANON}.gene_transcript_id.sorted.txt | join - ${CANON}.gene_type.txt | awk '{if ($3 == "protein_coding" || $3 == "pseudogene") print $2}' | sort \
        | join -1 1 -2 4 - <(sort -k 4,4 ${CANON}.extended_no_ovlp.3prime_extended.rnaseq_corrected.bed) \
        | awk 'BEGIN{OFS="\t";}{print $2, $3, $4, $1, $5, $6, $7, $8, $9, $10, $11, $12}' > ${CODING_PSEUDO_BED}
fi

[[ -s ${CANON}.extended_no_ovlp.3prime_extended.rnaseq_corrected.gtf ]] || die "extended GTF was not produced"
log "extended GTF: ${CANON}.extended_no_ovlp.3prime_extended.rnaseq_corrected.gtf ($(awk '$2 != 0' $ext/$base.filtered.3prime_extensions.txt | wc -l) transcripts extended)"
