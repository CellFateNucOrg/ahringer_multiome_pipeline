#!/bin/bash
# Stage 0: one-time setup. Run on the cluster LOGIN node (needs internet):
#   bash setup.sh            # everything
#   bash setup.sh envs       # only create conda environments
#   bash setup.sh project    # only create PROJECT_DIR layout, links and 10x whitelist
#   bash setup.sh reference  # only download the WormBase reference files
set -euo pipefail
PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PIPELINE_DIR
source "${PIPELINE_DIR}/scripts/common.sh"

what=${1:-all}

setup_envs() {
    command -v mamba >/dev/null || die "mamba not found in PATH"
    local key name
    for key in tools macs r meme; do
        case $key in
            tools) name=$ENV_TOOLS ;; macs) name=$ENV_MACS ;; r) name=$ENV_R ;; meme) name=$ENV_MEME ;;
        esac
        if conda env list | awk '{print $1}' | grep -qx "$name"; then
            log "conda env $name already exists"
        else
            log "creating conda env $name"
            mamba env create -y -n "$name" -f "${PIPELINE_DIR}/envs/mo_${key}.yaml"
        fi
    done
    activate_env "$ENV_R"
    Rscript "${PIPELINE_DIR}/envs/install_r_extras.R"
    # sinto from bioconda has a broken pkg_resources dependency; reinstall via pip to fix it
    activate_env "$ENV_TOOLS"
    pip install --force-reinstall --no-deps sinto --quiet
}

setup_project() {
    mkdir -p "$PROJECT_DIR"/{logs,plots,seurat_objects,star_idx,whitelist,cellranger_arc_library_files} \
             "$PROJECT_DIR"/data/sc/{cellranger_arc,rnaseq} "$PROJECT_DIR/$GENOME_DIR" "$PROJECT_DIR/$GA"
    ln -sfn "${PIPELINE_DIR}/external_data" "$PROJECT_DIR/data/external_data"
    ln -sfn "${PIPELINE_DIR}/scripts" "$PROJECT_DIR/scripts"

    [[ -x "${CELLRANGER_ARC_DIR}/cellranger-arc" ]] || die "cellranger-arc not found in CELLRANGER_ARC_DIR=${CELLRANGER_ARC_DIR}"
    # GEX barcode whitelist used by STARsolo (GEX barcodes are also the ones reported in the ATAC fragment file)
    local wl
    wl=$(find "$CELLRANGER_ARC_DIR" -path "*cellranger/barcodes/737K-arc-v1.txt.gz" 2>/dev/null | head -n 1)
    [[ -n "$wl" ]] || die "737K-arc-v1.txt.gz (GEX whitelist) not found under ${CELLRANGER_ARC_DIR}"
    zcat "$wl" > "$PROJECT_DIR/whitelist/737K-arc-v1.txt"
    log "whitelist: $wl -> $PROJECT_DIR/whitelist/737K-arc-v1.txt ($(wc -l < "$PROJECT_DIR/whitelist/737K-arc-v1.txt") barcodes)"
}

setup_reference() {
    cd "$PROJECT_DIR"
    local base="c_elegans.PRJNA13758.${WS}"
    wget -c -O "${CANON}.gtf.gz" "${WORMBASE_URL}/${base}.canonical_geneset.gtf.gz"
    wget -c -O "${GENOME_DIR}/elegans.fa.gz" "${WORMBASE_URL}/${base}.genomic.fa.gz"
    wget -c -O "${GA}/${base}.annotations.gff3.gz" "${WORMBASE_URL}/${base}.annotations.gff3.gz"
    gunzip -kf "${CANON}.gtf.gz"
    gunzip -kf "${GENOME_DIR}/elegans.fa.gz"
    log "reference files downloaded to $PROJECT_DIR"
}

case $what in
    envs) setup_envs ;;
    project) setup_project ;;
    reference) setup_reference ;;
    all) setup_envs; setup_project; setup_reference ;;
    *) die "unknown target '$what' (envs|project|reference|all)" ;;
esac
log "setup '$what' finished"
