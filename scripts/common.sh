#!/bin/bash
# Shared helpers, sourced by run_pipeline.sh and by every job script.
# Requires PIPELINE_DIR to be set.

[[ -n "${PIPELINE_DIR:-}" ]] || { echo "ERROR: PIPELINE_DIR is not set" >&2; exit 1; }
source "${PIPELINE_DIR}/config/config.sh"

die() { echo "ERROR: $*" >&2; exit 1; }
log() { echo "[$(date '+%F %T')] $*"; }

activate_env() {
    set +u
    if [[ -z "${CONDA_BASE:-}" ]]; then
        CONDA_BASE=$(conda info --base 2>/dev/null) || die "conda not found; set CONDA_BASE in config.sh"
    fi
    # shellcheck disable=SC1091
    source "${CONDA_BASE}/etc/profile.d/conda.sh"
    conda activate "$1" || die "cannot activate conda env $1"
    set -u
}

# ----- sample sheet helpers (tab-separated, header line, '#' comments allowed) -----
_samples_clean() { sed -e 's/\r$//' -e '/^[[:space:]]*#/d' -e '/^[[:space:]]*$/d' "$SAMPLES_TSV"; }

sample_ids() { _samples_clean | awk -F'\t' 'NR>1{print $1}'; }

n_samples() { sample_ids | wc -l; }

sample_by_index() { sample_ids | sed -n "${1}p"; }

# sample_get <sample_id> <column_name>
sample_get() {
    _samples_clean | awk -F'\t' -v id="$1" -v col="$2" '
        NR==1 { for (i=1;i<=NF;i++) h[$i]=i; if (!(col in h)) { print "missing column " col > "/dev/stderr"; exit 2 }; next }
        $1==id { print $(h[col]); found=1; exit }
        END { if (!found) exit 3 }'
}

# samples_where <column> <value> : sample ids whose column equals value
samples_where() {
    _samples_clean | awk -F'\t' -v col="$1" -v val="$2" '
        NR==1 { for (i=1;i<=NF;i++) h[$i]=i; next } $(h[col])==val { print $1 }'
}

# ----- derived paths (relative to PROJECT_DIR, identical to the paper's layout) -----
ANNOTATION=${WS}_extended_no_ovlp
GENOME_DIR=species/elegans/genome
GA=species/elegans/gene_annotation
CANON=${GA}/c_elegans.PRJNA13758.${WS}.canonical_geneset
CHR_SIZES=${GENOME_DIR}/elegans.chrom.sizes.txt
BLACKLIST=data/external_data/ce11-blacklist.bed
CR_REF=${WS}_cellranger_arc
BULK_PEAKS_DIR=scATAC_bulk_peaks
CODING_PSEUDO_BED=${CANON}.extended_no_ovlp.3prime_extended.rnaseq_corrected.coding_pseudogenes.bed
UTR5_BED=${CANON}.filtered.5prime_UTR_and_pre_coding_start.bed
ROUND1_DIR=cluster_specific_peaks/round_1_all_samples/${ANNOTATION}
ROUND2_DIR=cluster_specific_peaks/round_2_all_samples/${ANNOTATION}
MERGED_PREFIX=seurat_objects/combined_all_samples_SCT.${ANNOTATION}
FINAL_PREFIX=seurat_objects/wt_all.all_samples.all_cells.${ANNOTATION}

export WS ANNOTATION PROJECT_DIR PIPELINE_DIR
