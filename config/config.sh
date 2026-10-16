#!/bin/bash
# =====================================================================
#  Configuration for the C. elegans snMultiome pipeline
#  (reproduction of Ahringer lab scripts, bioRxiv 10.1101/2024.12.02.626321)
#  This file is sourced by every job. Edit the values below.
# =====================================================================

# ---------- Paths ----------
# Output/working directory on the cluster (created if missing). All results go here.
PROJECT_DIR=/mnt/meister.data/jsemple/20260728_ma_10x_multi_l3_PMW941_1/ahringer_pipeline_claude

# Sample sheet (tab-separated, see config/samples.tsv)
SAMPLES_TSV=${PIPELINE_DIR}/config/samples.tsv

# Cell Ranger ARC install root (the folder containing the 'cellranger-arc' executable).
# Download from https://www.10xgenomics.com/support/software/cell-ranger-arc/downloads
CELLRANGER_ARC_DIR=/mnt/meister.data/sharedSoftware/cellranger-arc-2.2.0

# ---------- Conda / mamba ----------
# Leave empty to auto-detect with 'conda info --base'
CONDA_BASE=""
ENV_TOOLS=mo_tools
ENV_MACS=mo_macs
ENV_R=mo_r
ENV_MEME=mo_meme

# ---------- Slurm ----------
#SLURM_ACCOUNT="jsemple"
#SLURM_PARTITION="all"
SLURM_EXTRA="--mail-type=FAIL --mail-user=jennifer.semple@unibe.ch"   # e.g. "--mail-type=END,FAIL --mail-user=you@unibe.ch"
SCRATCH_DIR="/scratch/${USER}"  # node-local scratch for BAM sorting temp files; leave empty to use the output dir

# Resources per job type (passed to sbatch; override freely)
RES_default="--cpus-per-task=1 --mem=8G --time=02:00:00"
RES_annotation="--cpus-per-task=2 --mem=8G --time=02:00:00"
RES_mkref="--cpus-per-task=12 --mem=64G --time=12:00:00"
RES_star_index="--cpus-per-task=12 --mem=48G --time=04:00:00"
RES_cellranger="--cpus-per-task=16 --mem=128G --time=48:00:00"
RES_fragments="--cpus-per-task=2 --mem=32G --time=06:00:00"
RES_starsolo="--cpus-per-task=12 --mem=64G --time=36:00:00"
RES_star_sort="--cpus-per-task=8 --mem=16G --time=04:00:00"
RES_dedup="--cpus-per-task=8 --mem=192G --time=48:00:00"
RES_bulk_peaks="--cpus-per-task=2 --mem=96G --time=12:00:00"
RES_utr_extension="--cpus-per-task=2 --mem=128G --time=12:00:00"
RES_r_small="--cpus-per-task=1 --mem=16G --time=04:00:00"
RES_emptydrops="--cpus-per-task=1 --mem=32G --time=06:00:00"
RES_soupx="--cpus-per-task=8 --mem=96G --time=12:00:00"
RES_merge="--cpus-per-task=4 --mem=192G --time=24:00:00"
RES_macs_group="--cpus-per-task=8 --mem=32G --time=06:00:00"
RES_bed_merge="--cpus-per-task=1 --mem=16G --time=04:00:00"
RES_r_large="--cpus-per-task=4 --mem=192G --time=48:00:00"
RES_idr="--cpus-per-task=8 --mem=32G --time=04:00:00"
RES_rna_tracks="--cpus-per-task=12 --mem=64G --time=24:00:00"
RES_meme="--cpus-per-task=8 --mem=32G --time=24:00:00"
# diagnostics (neuron_diag, neuron_atac_diag, atac_qc_diag). Measured for 1 sample (14,265 nuclei):
# atac_qc_diag peaked at 9.5 GB and took 58 min. Scale memory and time roughly with the number of samples
# (1 sample: --mem=16G --time=02:00:00 is enough; 4 samples: the values below).
RES_diag="--cpus-per-task=2 --mem=32G --time=04:00:00"

# ---------- Reference ----------
WS=WS285                  # WormBase release used in the paper
# EBI mirror of WormBase (downloads.wormbase.org may block scripted downloads)
WORMBASE_URL=https://ftp.ebi.ac.uk/pub/databases/wormbase/releases/${WS}/species/c_elegans/PRJNA13758

# ---------- Parameters (paper defaults) ----------
# per-sample QC (seurat_soupx)
MIN_ATAC_FRAGS=200        # min ATAC counts in peaks per barcode
MT_PCT_MAX=5              # max % mitochondrial UMIs
MAX_RNA_UMI=10000         # max RNA UMIs per nucleus
SOUPX_QUANTILE=0.80
# L3 gene sets for per-sample QC (seurat_soupx.R)
L3_MARKER_RANKED_FILE=${PIPELINE_DIR}/tissueMarkers/consensus/markers_consensus_min2sources.csv
GERMLINE_TISSUE=germline        # tissue of L3_MARKER_RANKED_FILE whose clusters are excluded from the SoupX estimate
N_MARKERS_PER_TISSUE_PLOT=3     # top markers per tissue shown before / after SoupX
N_SOUP_CHANGE_GENES=6           # tissue-specific soup genes shown in the SoupX change maps
N_SOUP_SCORE_GENES=50           # top soup genes averaged into the per-nucleus soup score

# STARsolo: extra options, e.g. "--soloBarcodeReadLength 0" if GEX R1 is longer than 28 bp
STAR_EXTRA_OPTS="--soloBarcodeReadLength 0"
DEDUP_PARALLEL=2          # reduced from 4: each split BAM ~3.4GB, umi_tools uses ~10x RAM

# 3' UTR extension
UTR_EXT_BIN=100
UTR_EXT_DISTANCES=0,100,200,300,400

# batch-enriched cluster removal (batch_enriched_clusters_removal.R)
BATCH_FOLD=4                        # flag cluster if one batch > FOLD x another batch of the same stage group (share of each batch's own nuclei)
LATE_STAGE_GROUPS="late"            # stage groups for which clusters are flagged when they contribute > LATE_MAX_FRACTION
LATE_MAX_FRACTION=0.8
MIN_WNN_FEATURES=200
BATCH_MIN_GROUP_NUCLEI=30           # a stage group is only compared for a cluster if it has >= this many nuclei in it ...
BATCH_MIN_GROUP_FRACTION=0.2        # ... and makes up >= this fraction of the cluster
BATCH_REMOVAL_MODE=remove           # remove = drop flagged clusters; report = only flag them (table + plots), keep all nuclei
# genes whose clusters are never removed as batch-enriched (paper: cpg-3 ccch-2 for 2-cell / ABx embryo cells);
# space-separated, empty = no protection
PROTECT_CLUSTER_MARKERS=""

# final cell type assignment
# optional table "cluster_number<TAB>lineage_label" for late clusters (see config/late_cluster_reassignment.tsv)
LATE_CLUSTER_TABLE=${PIPELINE_DIR}/config/late_cluster_reassignment.tsv
# marker tables for L3 tissue annotation (l3_lineage_annotation.R); one TSV path per line
MARKER_LIST_FILE=${PIPELINE_DIR}/tissueMarkers/pipeline_markers_v2/marker_list.txt
# optional taxonomy TSV for raw->tissue_type lookup (for per-study marker files)
MARKER_TAXONOMY_FILE=${PIPELINE_DIR}/cellMarkers/cell_type_taxonomy.tsv
# tables with gene_label + WBGeneID columns, used to translate marker labels to ${WS} gene names
MARKER_GENE_ID_TABLES="${PIPELINE_DIR}/tissueMarkers/evidence/marker_evidence_long_tissue.csv ${PIPELINE_DIR}/tissueMarkers/consensus/markers_consensus_min2sources.csv"

# cell types relabelled as 'unassigned' after inspection (paper: "Cxx"); space-separated, may be empty
UNASSIGN_CELL_TYPES="Cxx"

# IDR on pseudo-replicates (samples split with the 'idr_rep' column of samples.tsv)
IDR_FINAL_SET=0.05        # which IDR set builds the final peak set: 0.05 (paper) or 0.01
IDR_QC_HEATMAPS=yes

# optional motif analysis: MEME motif database (from https://meme-suite.org/meme/doc/download.html)
MEME_DB=${PROJECT_DIR}/motif_databases/CIS-BP_2.00/Caenorhabditis_elegans.meme
