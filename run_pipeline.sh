#!/bin/bash
# Submit one stage of the snMultiome pipeline to Slurm.
# Usage: ./run_pipeline.sh <stage> [--after JOBID[:JOBID...]] [--dry-run]
set -euo pipefail
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PIPELINE_DIR
source "${PIPELINE_DIR}/scripts/common.sh"
J=${PIPELINE_DIR}/stages/jobs

usage() {
    cat <<EOF
Usage: $(basename "$0") <stage> [--after JOBID[:JOBID]] [--dry-run]

  1  reference        annotation derivatives, cellranger-arc mkref, STAR index (${WS})
  2  align            cellranger-arc count, fragment files, STARsolo (${WS}), UMI dedup
  3  utr_extension    bulk scATAC peaks, 3' UTR extension, STAR index (${ANNOTATION})
  4  starsolo         STARsolo on extended annotation + barcode-rank QC plots
                      >>> MANUAL: set 'barhop' and 'ambient' per sample in samples.tsv
  5  sample_qc        EmptyDrops, per-sample Seurat / scDblFinder / SoupX / ATAC
  6  merge            merge samples, RNA clustering, per-cluster fragments
  7  wnn_annotation   MACS2 per cluster, WNN, batch-enriched cluster removal, marker-based annotation
                      >>> MANUAL (optional): fill config/late_cluster_reassignment.tsv
  8  cell_types       final cell type / lineage assignment
  9  peaks_idr        MACS2 per cell type & lineage, IDR, final peak set, FRiP, cell-cycle classifier
  10 tracks           (optional) RNA and ATAC bigWig tracks per cell type / embryo stage
  11 transient_genes  (optional) pseudo-bulk DESeq2 of pre- vs post-gastrulation lineages
  motifs <bed> <name> (optional) MEME-ChIP + FIMO on any BED file

Options:
  --after IDS   make every job of this stage depend on these job ids (afterok)
  --dry-run     print the sbatch commands without submitting
EOF
    exit 1
}

[[ $# -ge 1 ]] || usage
STAGE=$1; shift
AFTER=""; DRY=0; POS=()
while [[ $# -gt 0 ]]; do
    case $1 in
        --after) AFTER=$2; shift 2 ;;
        --dry-run) DRY=1; shift ;;
        -h|--help) usage ;;
        *) POS+=("$1"); shift ;;
    esac
done

[[ -d "$PROJECT_DIR" ]] || die "PROJECT_DIR $PROJECT_DIR does not exist; run 'bash setup.sh' first"
mkdir -p "$PROJECT_DIR/logs"
cd "$PROJECT_DIR"

# submit <res_key> <name> <deps or ""> [extra sbatch options] -- <script> [args]
# prints the job id on stdout
submit() {
    local res_key=$1 name=$2 dep=$3; shift 3
    local -a opts=()
    while [[ $# -gt 0 && $1 != "--" ]]; do opts+=("$1"); shift; done
    [[ ${1:-} == "--" ]] && shift
    local res_var="RES_${res_key}" log_pat="%x.%j.out"
    local -a res extra
    read -r -a res <<< "${!res_var:-$RES_default}"
    read -r -a extra <<< "${SLURM_EXTRA}"
    [[ " ${opts[*]:-} " == *" --array="* ]] && log_pat="%x.%A_%a.out"
    local deps=${dep}
    [[ -n "$AFTER" ]] && deps=${deps:+${deps}:}${AFTER}
    local -a cmd=(sbatch --parsable --job-name="mo_${name}" --output="$PROJECT_DIR/logs/${log_pat}"
                  --export="ALL,PIPELINE_DIR=${PIPELINE_DIR}")
    [[ -n "$deps" ]] && cmd+=(--dependency="afterok:${deps}" --kill-on-invalid-dep=yes)
    cmd+=("${res[@]}" "${extra[@]}" "${opts[@]}" "$@")
    if [[ $DRY -eq 1 ]]; then
        echo "[dry-run] ${cmd[*]}" >&2
        echo "DRY_${name}"
    else
        local id
        id=$("${cmd[@]}" | cut -d';' -f1)
        echo "submitted mo_${name}: job ${id}${deps:+ (after ${deps})}" >&2
        echo "$id"
    fi
}

need() { # need <path> <message>
    if [[ ! -e "$1" ]]; then
        if [[ $DRY -eq 1 ]]; then echo "[dry-run] WARNING missing: $1 ($2)" >&2; else die "missing $1 -- $2"; fi
    fi
}

n_groups() { # number of groups in a cells_per_* table (header line + one line per group)
    if [[ -f "$1" ]]; then echo $(( $(wc -l < "$1") - 1 )); else echo 1; fi
}

N=$(n_samples)
[[ $N -ge 1 ]] || die "no samples in $SAMPLES_TSV"

case $STAGE in
1|reference)
    need "${CANON}.gtf" "run 'bash setup.sh reference'"
    need "${GENOME_DIR}/elegans.fa" "run 'bash setup.sh reference'"
    a=$(submit annotation ref_annotation "" -- "$J/ref_annotation.sh")
    submit mkref cr_mkref "" -- "$J/cellranger_mkref.sh" >/dev/null
    submit star_index star_index_${WS} "$a" -- "$J/star_index.sh" "${CANON}.gtf" "$WS" >/dev/null
    ;;
2|align)
    need "whitelist/737K-arc-v1.txt" "run 'bash setup.sh project'"
    need "$CR_REF/fasta/genome.fa" "stage 1 (cellranger-arc mkref)"
    need "star_idx/$WS/SA" "stage 1 (STAR index)"
    c=$(submit cellranger cr_count "" --array=1-"$N" -- "$J/cellranger_count.sh")
    submit fragments atac_fragments "$c" --array=1-"$N" -- "$J/atac_fragments.sh" >/dev/null
    s=$(submit starsolo starsolo_${WS} "" --array=1-"$N" -- "$J/starsolo.sh" "$WS")
    submit dedup star_dedup "$s" --array=1-"$N" -- "$J/star_dedup.sh" "$WS" >/dev/null
    ;;
3|utr_extension)
    [[ -n "$(samples_where utr_ext yes)" ]] || die "no sample has utr_ext=yes in samples.tsv"
    for id in $(samples_where utr_ext yes); do
        need "data/sc/cellranger_arc/$id/outs/atac_fragments.tsv.gz" "stage 2"
        need "data/sc/rnaseq/$id/star_${WS}/star_${WS}_Aligned.sortedByCoord.out.split.merge.dedup.bam" "stage 2"
    done
    b=$(submit bulk_peaks bulk_peaks "" -- "$J/bulk_peaks.sh")
    u=$(submit utr_extension utr_extension "$b" -- "$J/utr_extension.sh")
    submit star_index star_index_ext "$u" -- "$J/star_index.sh" \
        "${CANON}.extended_no_ovlp.3prime_extended.rnaseq_corrected.gtf" "$ANNOTATION" >/dev/null
    ;;
4|starsolo)
    need "star_idx/$ANNOTATION/SA" "stage 3"
    s=$(submit starsolo starsolo_ext "" --array=1-"$N" -- "$J/starsolo.sh" "$ANNOTATION")
    submit r_small barcode_qc "$s" --array=1-"$N" -- "$J/r_script.sh" barcode_rank_qc.R __SAMPLE__ "$ANNOTATION" >/dev/null
    echo "When finished: inspect $PROJECT_DIR/plots/barcode_rank_qc/*.pdf and set 'barhop' and 'ambient' in $SAMPLES_TSV" >&2
    ;;
5|sample_qc)
    for id in $(sample_ids); do
        for col in barhop ambient; do
            v=$(sample_get "$id" "$col")
            [[ $v =~ ^[0-9]+$ ]] || die "sample $id: '$col' must be an integer (is '$v'); see stage 4 QC plots"
        done
        need "data/sc/rnaseq/$id/star_${ANNOTATION}/star_${ANNOTATION}_Solo.out/GeneFull/raw/um_reads/matrix.mtx.gz" "stage 4"
        need "data/sc/cellranger_arc/$id/outs/atac_fragments_no_dash.tsv.gz" "stage 2"
    done
    need "$BULK_PEAKS_DIR/macs/resized_fragments.all_samples_peaks.bed" "stage 3"
    e=$(submit emptydrops emptydrops "" --array=1-"$N" -- "$J/r_script.sh" emptydrops.R __SAMPLE__ "$ANNOTATION")
    submit soupx seurat_soupx "$e" --array=1-"$N" -- "$J/r_script.sh" seurat_soupx.R __SAMPLE__ "$ANNOTATION" \
        "$MIN_ATAC_FRAGS" "$BULK_PEAKS_DIR/macs/resized_fragments.all_samples_peaks.bed" >/dev/null
    ;;
6|merge)
    for id in $(sample_ids); do
        need "seurat_objects/$id.$ANNOTATION.postSoupX.ATAC_MACS.rds" "stage 5"
        need "data/sc/cellranger_arc/$id/outs/atac_fragments_barcode_corrected.tsv.gz" "stage 2"
    done
    submit merge merge_samples "" -- "$J/r_script.sh" merge_multiple_seurat_objects.R \
        "$(sample_ids | paste -sd,)" "$BULK_PEAKS_DIR/macs/resized_fragments.all_samples_peaks.bed" \
        "$ROUND1_DIR/bed" "$ANNOTATION" >/dev/null
    ;;
7|wnn_annotation)
    cpc="${MERGED_PREFIX}.all_cells.MACS_peaks.cells_per_cluster.txt"
    need "$cpc" "stage 6"
    m=$(submit macs_group macs_round1 "" --array=1-"$(n_groups "$cpc")" -- "$J/macs_group.sh" "$ROUND1_DIR" "$cpc")
    p=$(submit bed_merge peaks_round1 "$m" -- "$J/merge_group_peaks.sh" round1)
    w=$(submit r_large wnn "$p" -- "$J/r_script.sh" merge_multiple_seurat_objects.cluster_peaks_step_2.R \
        "${MERGED_PREFIX}.all_cells.MACS_peaks.rds" "$ROUND1_DIR/all_peaks_merge.intergenic.00001.bed" WNN "$ANNOTATION")
    b=$(submit r_large batch_removal "$w" -- "$J/r_script.sh" batch_enriched_clusters_removal.R \
        "${MERGED_PREFIX}.WNN.rds" "$ANNOTATION")
    submit r_large lineage_annotation "$b" -- "$J/r_script.sh" lineage_specific_annotation.R \
        "${MERGED_PREFIX}.WNN_clean.rds" cell_type_annotation_all "$ANNOTATION" >/dev/null
    echo "When finished: inspect plots/lineage_annotation and optionally fill $LATE_CLUSTER_TABLE" >&2
    ;;
8|cell_types)
    need "${MERGED_PREFIX}.WNN_clean.rds" "stage 7"
    need "cell_type_annotation_all/AB.cell_assignment.$ANNOTATION.txt" "stage 7"
    submit r_large final_cell_types "" -- "$J/r_script.sh" final_cell_type_assignment.R \
        "${MERGED_PREFIX}.WNN_clean.rds" "$ANNOTATION" "$LATE_CLUSTER_TABLE" >/dev/null
    ;;
9|peaks_idr)
    cpl="${FINAL_PREFIX}.cells_per_cell_lineage.txt"; cpt="${FINAL_PREFIX}.cells_per_cell_type.txt"
    need "$cpl" "stage 8"; need "$cpt" "stage 8"
    ml=$(submit macs_group macs_lineage "" --array=1-"$(n_groups "$cpl")" -- "$J/macs_group.sh" "$ROUND2_DIR/cell_lineage" "$cpl")
    mt=$(submit macs_group macs_cell_type "" --array=1-"$(n_groups "$cpt")" -- "$J/macs_group.sh" "$ROUND2_DIR/cell_type" "$cpt")
    p=$(submit bed_merge peaks_round2 "${ml}:${mt}" -- "$J/merge_group_peaks.sh" round2)
    f=$(submit r_large frags_for_idr "$p" -- "$J/r_script.sh" ATAC_frags_for_IDR.R \
        "${FINAL_PREFIX}.rds" "$ROUND2_DIR/all_peaks.pre_IDR.bed" "$ROUND2_DIR/cell_type/macs_peak_frag_counts")
    d=$(submit idr idr "$f" --array=1-"$(n_groups "$cpt")" -- "$J/idr_atac.sh" "$cpt")
    fp=$(submit bed_merge final_peakset "$d" -- "$J/final_peakset.sh")
    ip=$(submit r_large integrate_peaks "$fp" -- "$J/r_script.sh" seurat_integrate_cluster_peaks.R \
        "${FINAL_PREFIX}.rds" "$ROUND2_DIR/all_peaks_merge.intergenic.00001.bed" final_peakset \
        "${FINAL_PREFIX}.final_peakset_IDR.rds" "$ANNOTATION")
    submit r_large cell_cycle "$ip" -- "$J/r_script.sh" cell_cycle_random_forest.R \
        "${FINAL_PREFIX}.final_peakset_IDR.rds" "${FINAL_PREFIX}.final_peakset_IDR_cell_cycle.rds" >/dev/null
    ;;
10|tracks)
    cpt="${FINAL_PREFIX}.cells_per_cell_type.txt"
    need "${FINAL_PREFIX}.final_peakset_IDR.rds" "stage 9"
    t=$(submit r_small track_tables "" -- "$J/r_script.sh" export_track_tables.R "${FINAL_PREFIX}.final_peakset_IDR.rds" tracks)
    s=$(submit rna_tracks rna_split "$t" --array=1-"$N" -- "$J/rna_split_by_cell_type.sh")
    submit rna_tracks rna_cell_type_bw "$s" --array=1-"$(n_groups "$cpt")" -- "$J/rna_cell_type_tracks.sh" "$cpt" >/dev/null
    submit rna_tracks atac_stage_bw "$t" -- "$J/atac_stage_tracks.sh" >/dev/null
    ;;
11|transient_genes)
    need "${FINAL_PREFIX}.final_peakset_IDR_cell_cycle.rds" "stage 9"
    submit r_large transient_genes "" -- "$J/r_script.sh" transient_genes_DESeq.R \
        "${FINAL_PREFIX}.final_peakset_IDR_cell_cycle.rds" >/dev/null
    ;;
motifs)
    [[ ${#POS[@]} -eq 2 ]] || die "usage: $(basename "$0") motifs <bed file> <output name>"
    need "${POS[0]}" "input BED"
    need "$ROUND2_DIR/all_peaks.bed" "stage 9 (background model)"
    submit meme motifs_"${POS[1]}" "" -- "$J/motifs.sh" "$(readlink -f "${POS[0]}")" "${POS[1]}" >/dev/null
    ;;
*) usage ;;
esac
