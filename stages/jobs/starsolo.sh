#!/bin/bash
# STARsolo on one sample (STAR_launcher.sh + STARsolo_output_parser). Args: <index name>
set -euo pipefail
source "${PIPELINE_DIR}/scripts/common.sh"
cd "$PROJECT_DIR"
activate_env "$ENV_TOOLS"

idx=$1
id=$(sample_by_index "${SLURM_ARRAY_TASK_ID}")
outdir=data/sc/rnaseq/$id/star_$idx
prefix=$outdir/star_${idx}_
mkdir -p $outdir

# ── Skip if complete ─────────────────────────────────────────────────────────
if [[ -s ${prefix}Solo.out/GeneFull/raw/UniqueAndMult-EM.mtx ]]; then
    log "$id: STARsolo output exists, skipping alignment"
else
    dir=$(sample_get "$id" gex_fastq_dir)
    r1=(); r2=()
    IFS=, read -r -a prefixes <<< "$(sample_get "$id" gex_fastq_prefix)"
    for p in "${prefixes[@]}"; do
        mapfile -t f1 < <(ls "$dir"/${p}_S*_R1_001.fastq.gz | sort)
        mapfile -t f2 < <(ls "$dir"/${p}_S*_R2_001.fastq.gz | sort)
        [[ ${#f1[@]} -gt 0 && ${#f1[@]} -eq ${#f2[@]} ]] || die "$id: R1/R2 FASTQs not found or unpaired for prefix $p in $dir"
        r1+=("${f1[@]}"); r2+=("${f2[@]}")
    done
    R1=$(IFS=,; echo "${r1[*]}"); R2=$(IFS=,; echo "${r2[*]}")
    log "$id: R2=$R2"
    log "$id: R1=$R1"

    # ── Use local scratch if available (faster I/O, clean failure recovery) ──
    # Scratch dir is named after sample+index (not job ID) so every re-run of the same
    # sample wipes any stale dir left behind by a SIGKILL'd previous attempt.
    if [[ -n "${SCRATCH_DIR:-}" ]]; then
        scratch="${SCRATCH_DIR}/starsolo_${id}_${idx}"
        # Wipe any stale scratch from a previously SIGKILL'd run
        if [[ -d "$scratch" ]]; then
            log "$id: removing stale scratch dir $scratch from a previous (killed?) run"
            rm -rf "$scratch"
        fi
        mkdir -p "$scratch"
        run_prefix="$scratch/star_${idx}_"
        # On graceful failure (set -e / SIGTERM): remove scratch dir
        # SIGKILL cannot be trapped; the fresh rm -rf above handles that on next run
        trap 'log "$id: STAR failed or interrupted; cleaning scratch $scratch"; rm -rf "$scratch"' ERR EXIT
        log "$id: writing STAR output to scratch $scratch"
    else
        run_prefix="$prefix"
        # On failure: remove any partial STAR output from NFS so next run starts clean
        trap 'log "$id: STAR failed or interrupted; removing partial output under $outdir"; \
              rm -rf "${prefix}Solo.out" "${prefix}_STARtmp" "${prefix}Aligned.out.bam"' ERR EXIT
        log "$id: no scratch dir; writing STAR output to NFS directly"
    fi

    # shellcheck disable=SC2086
    STAR --genomeDir star_idx/$idx --readFilesIn $R2 $R1 --soloType CB_UMI_Simple \
         --soloCBwhitelist whitelist/737K-arc-v1.txt --soloUMIlen 12 --soloCellFilter None --soloFeatures GeneFull \
         --soloMultiMappers EM --readFilesCommand zcat --runThreadN "${SLURM_CPUS_PER_TASK:-12}" \
         --outFileNamePrefix "$run_prefix" --outSAMtype BAM Unsorted --twopassMode None \
         --outSAMattributes NH HI nM AS CR UR GX GN sS sQ sM RG XS \
         --outSAMattrRGline "ID:"$id --alignIntronMax 20000 ${STAR_EXTRA_OPTS}

    # ── Copy scratch → NFS and disarm the failure trap ────────────────────────
    if [[ -n "${SCRATCH_DIR:-}" ]]; then
        trap - ERR EXIT  # disarm before copy; copy failure leaves scratch intact for inspection
        log "$id: copying STAR outputs from scratch to $outdir"
        cp "${run_prefix}Aligned.out.bam" "${prefix}Aligned.out.bam"
        cp -r "${run_prefix}Solo.out" "${prefix}Solo.out"
        for f in Log.out Log.final.out Log.progress.out SJ.out.tab; do
            [[ -e "${run_prefix}${f}" ]] && cp "${run_prefix}${f}" "${prefix}${f}"
        done
        rm -rf "$scratch"
        log "$id: scratch copy complete; scratch removed"
    else
        trap - ERR EXIT  # STAR completed; disarm cleanup trap
    fi
fi

# ── Parser ───────────────────────────────────────────────────────────────────
if [[ -s ${prefix}Solo.out/GeneFull/raw/um_reads/matrix.mtx.gz ]]; then
    log "$id: STARsolo parser output exists, skipping"
else
    python scripts/STARsolo_output_parser.py $id $idx ${CANON}.gene_name.txt
fi
log "$id: STARsolo ($idx) done"
