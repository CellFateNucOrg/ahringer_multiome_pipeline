#!/bin/bash
# Gene/genome annotation preprocessing (sc_project_pipeline.20241202.txt, "gene/genome annotation preprocessing")
set -euo pipefail
source "${PIPELINE_DIR}/scripts/common.sh"
cd "$PROJECT_DIR"
activate_env "$ENV_TOOLS"

log "genome index and chromosome sizes"
samtools faidx ${GENOME_DIR}/elegans.fa
cut -f 1,2 ${GENOME_DIR}/elegans.fa.fai > ${CHR_SIZES}
sort -k 1,1 ${CHR_SIZES} > ${GENOME_DIR}/elegans.chrom.sizes.sorted.txt
awk 'BEGIN{OFS="\t";}{print $1, 1, $2}' ${CHR_SIZES} > ${GENOME_DIR}/elegans.bed

log "GTF -> BED12"
gtfToGenePred ${CANON}.gtf ${CANON}.genePred
genePredToBed ${CANON}.genePred ${CANON}.bed.temp
sort -k 4,4 ${CANON}.bed.temp > ${CANON}.bed
rm ${CANON}.bed.temp ${CANON}.genePred

log "gene type / name tables"
awk 'BEGIN{FS="\t|\""; OFS="\t";}{if ($3 == "gene") print $10, $16}' ${CANON}.gtf | sort -k 1,1 > ${CANON}.gene_type.txt
awk 'BEGIN{FS="\t|\""; OFS="\t";}{if ($3 == "gene") print $18, $16}' ${CANON}.gtf | sort -k 1,1 > ${CANON}.gene_name_gene_type.txt

python scripts/gtf_gene_type_filter.py ${CANON}.gtf lincRNA,protein_coding,pseudogene ${CANON}.filtered.gtf ${CANON}.gene_type.txt

awk 'BEGIN{FS="\t|\""; OFS="\t";}{if ($3 == "gene") print $10, $18}' ${CANON}.filtered.gtf | sort -k 1,1 > ${CANON}.filtered.gene_name.txt
awk 'BEGIN{FS="\t|\""; OFS="\t";}{if ($3 == "gene") print $10, $18}' ${CANON}.gtf | sort -k 1,1 > ${CANON}.gene_name.txt

# gene to transcript tables
awk 'BEGIN{FS="\t|\""; OFS="\t";}{if ($3 == "transcript") print $10, $14}' ${CANON}.filtered.gtf | sort -k 2,2 > ${CANON}.filtered.gene_transcript_id.sorted.txt
awk 'BEGIN{FS="\t|\""; OFS="\t";}{if ($3 == "transcript") print $10, $14}' ${CANON}.gtf | sort -k 2,2 > ${CANON}.gene_transcript_id.sorted.txt
sort -k 1,1 ${CANON}.gene_transcript_id.sorted.txt > ${CANON}.gene_transcript_id.WB_sorted.txt

# MtDNA genes (common names)
awk 'BEGIN{FS="\t|\""; OFS="\t";}{if ($3 == "gene" && $1 == "MtDNA") print $10}' ${CANON}.gtf | sort | join - ${CANON}.filtered.gene_name.txt | awk '{print $2}' | sort > ${CANON}.filtered.MtDNA_genes.genes

# BED12 of the filtered gene set (the paper creates it inside gtf_3_utr_extension.20231108.py; identical procedure)
gtfToGenePred ${CANON}.filtered.gtf ${CANON}.filtered.genePred
genePredToBed ${CANON}.filtered.genePred ${CANON}.filtered.bed.temp
awk 'BEGIN{FS="\t";OFS="\t";}{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' ${CANON}.filtered.bed.temp | sort -k 4,4 > ${CANON}.filtered.bed
rm ${CANON}.filtered.bed.temp

# putative promoter region (-1000/+200bp from TSS) of filtered genes
flankBed -l 1000 -r 0 -s -i ${CANON}.filtered.bed -g ${CHR_SIZES} | cut -f 1,2,3,4,5,6 | slopBed -l 0 -r 200 -s -i stdin -g ${CHR_SIZES} | sort -k 1,1 -k2,2n > ${CANON}.filtered.promoter_region.bed

# 5' UTRs of filtered genes (>=100bp) - used to filter accessible sites
awk 'BEGIN{OFS="\t";}{if ($6 == "+" && $7 == $8) print $1, $2, $2 + 100, $4, $5, $6; else if ($6 == "-" && $7 == $8) print $1, $3-100, $3, $4, $5, $6; else if ($6 == "+" && $7 - $2 >= 100) print $1, $2, $7, $4, $5, $6; else if ($6 == "-" && $3 - $8 >= 100) print $1, $8, $3, $4, $5, $6; else if ($6 == "+" && $7 - $2 < 100) print $1, $7 - 100, $7, $4, $5, $6; else if ($6 == "-" && $3 - $8 < 100) print $1, $8, $8 + 100, $4, $5, $6}' ${CANON}.filtered.bed > ${UTR5_BED}

# operons
gff3=${GA}/c_elegans.PRJNA13758.${WS}.annotations.gff3.gz
if [[ -s $gff3 ]]; then
    zcat $gff3 | grep operon | awk 'BEGIN{FS="\t|=|;"; OFS="\t";}{if ($2 == "dicistronic_mRNA" || $2 == "operon") print $1, $4, $5, $10, $12, $7}' > ${GA}/c_elegans.PRJNA13758.${WS}.operons.bed
    python scripts/operon_parser.py ${GA}/c_elegans.PRJNA13758.${WS}.operons.bed ${CANON}.bed ${CANON}.gene_transcript_id.sorted.txt ${GA}/c_elegans.PRJNA13758.${WS}.operons.gene_position.txt > ${GA}/operon_parser.log
fi

log "annotation done"
