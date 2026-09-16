#!/bin/bash
# cellranger-arc mkref on WormBase genome + canonical gene set
set -euo pipefail
source "${PIPELINE_DIR}/scripts/common.sh"
cd "$PROJECT_DIR"

if [[ -s "${CR_REF}/fasta/genome.fa" ]]; then
    log "${CR_REF} already exists, skipping"; exit 0
fi

cat > cellranger_arc_mkref.config <<EOF
{
    organism: "Caenorhabditis_elegans"
    genome: ["${CR_REF}"]
    input_fasta: ["${PROJECT_DIR}/${GENOME_DIR}/elegans.fa"]
    input_gtf: ["${PROJECT_DIR}/${CANON}.gtf"]
    non_nuclear_contigs: ["MtDNA"]
}
EOF

"${CELLRANGER_ARC_DIR}/cellranger-arc" mkref --config=cellranger_arc_mkref.config --nthreads="${SLURM_CPUS_PER_TASK:-12}"
log "mkref done"
