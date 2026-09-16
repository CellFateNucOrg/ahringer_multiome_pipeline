#!/bin/bash
# De novo motif discovery (memechip.20240521.sh) + logos + FIMO scan of all accessible sites. Args: <bed> <name>
set -euo pipefail
source "${PIPELINE_DIR}/scripts/common.sh"
cd "$PROJECT_DIR"
activate_env "$ENV_MEME"

input_bed=$1
name=$2
out=motifs/$name
mkdir -p $out
genome_fa=${GENOME_DIR}/elegans.fa

# 2nd-order Markov background from the final accessible-site set
[[ -s $ROUND2_DIR/all_peaks.fa ]] || fastaFromBed -bed $ROUND2_DIR/all_peaks.bed -fi $genome_fa -fo $ROUND2_DIR/all_peaks.fa
[[ -s $ROUND2_DIR/all_peaks.bkg ]] || fasta-get-markov -m 2 $ROUND2_DIR/all_peaks.fa $ROUND2_DIR/all_peaks.bkg

cut -f 1,2,3 $input_bed | sort -k 1,1 -k2,2n | uniq > $out/$name.loci.bed
fastaFromBed -fi $genome_fa -bed $out/$name.loci.bed -fo $out/$name.fa

min_w=6; max_w=15; n_motifs=20
db_opt=()
if [[ -s "$MEME_DB" ]]; then db_opt=(-db "$MEME_DB"); else log "WARNING: MEME_DB not found ($MEME_DB); running without a motif database"; fi
meme-chip -ccut 0 -oc $out/$name.w_${min_w}_${max_w} -minw $min_w -maxw $max_w -seed 32 -filter-thresh 0.05 -meme-mod zoops \
    -meme-nmotifs $n_motifs -meme-minsites 10 -meme-p "${SLURM_CPUS_PER_TASK:-8}" -streme-pvt 0.05 "${db_opt[@]}" \
    -bfile $ROUND2_DIR/all_peaks.bkg $out/$name.fa

res=$out/$name.w_${min_w}_${max_w}
grep MOTIF $res/combined.meme | awk 'BEGIN{FS=" |-"; OFS="\t";}{if ($4 == "MEME") print NR, $4 "-" $5; else print NR, $5 "-" $6}' > $res/combined.meme.motif_names
mkdir -p $res/logos
while read -r motif_n motif_name; do
    ceqlogo -i $res/combined.meme -m $motif_n -o $res/logos/$motif_name.fwd.png -f PNG
    ceqlogo -i $res/combined.meme -m $motif_n -o $res/logos/$motif_name.rev.png -f PNG -r
done < $res/combined.meme.motif_names

fimo --max-stored-scores 1000000 --oc $out/fimo --thresh 0.001 $res/combined.meme $ROUND2_DIR/all_peaks.fa
rm $out/$name.fa
log "motifs for $name in $out"
