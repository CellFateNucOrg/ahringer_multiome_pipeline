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

if [[ -s ${prefix}Solo.out/GeneFull/raw/UniqueAndMult-EM.mtx && -s ${prefix}Aligned.sortedByCoord.out.bam ]]; then
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
    # Use 75% of the job's allocated memory for BAM sorting (in bytes)
    bam_sort_ram=$(( ${SLURM_MEM_PER_NODE:-96000} * 1024 * 1024 * 3 / 4 ))
    # shellcheck disable=SC2086
    STAR --genomeDir star_idx/$idx --readFilesIn $R2 $R1 --soloType CB_UMI_Simple \
         --soloCBwhitelist whitelist/737K-arc-v1.txt --soloUMIlen 12 --soloCellFilter None --soloFeatures GeneFull \
         --soloMultiMappers EM --readFilesCommand zcat --runThreadN "${SLURM_CPUS_PER_TASK:-12}" \
         --outFileNamePrefix $prefix --outSAMtype BAM SortedByCoordinate --outWigType wiggle --twopassMode None \
         --limitBAMsortRAM $bam_sort_ram --outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM RG XS \
         --outSAMattrRGline "ID:"$id --alignIntronMax 20000 ${STAR_EXTRA_OPTS}
fi

python scripts/STARsolo_output_parser.py $id $idx ${CANON}.gene_name.txt
log "$id: STARsolo ($idx) done"
