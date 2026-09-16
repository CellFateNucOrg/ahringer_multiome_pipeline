#!/bin/bash
# =====================================================================
#  Configuration for the C. elegans snMultiome pipeline
#  (reproduction of Ahringer lab scripts, bioRxiv 10.1101/2024.12.02.626321)
#  This file is sourced by every job. Edit the values below.
# =====================================================================

# ---------- Paths ----------
# Output/working directory on the cluster (created if missing). All results go here.
PROJECT_DIR=/mnt/meister.data/jsemple/20260728_ma_10x_multi_PMW941_1/ahringer_pipeline_claude

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
SLURM_EXTRA=""            # e.g. "--mail-type=END,FAIL --mail-user=you@unibe.ch"

# Resources per job type (passed to sbatch; override freely)
RES_default="--cpus-per-task=1 --mem=8G --time=02:00:00"
RES_annotation="--cpus-per-task=2 --mem=8G --time=02:00:00"
RES_mkref="--cpus-per-task=12 --mem=64G --time=12:00:00"
RES_star_index="--cpus-per-task=12 --mem=48G --time=04:00:00"
RES_cellranger="--cpus-per-task=16 --mem=128G --time=48:00:00"
RES_fragments="--cpus-per-task=2 --mem=32G --time=06:00:00"
RES_starsolo="--cpus-per-task=12 --mem=96G --time=24:00:00"
RES_dedup="--cpus-per-task=8 --mem=96G --time=24:00:00"
RES_bulk_peaks="--cpus-per-task=2 --mem=48G --time=12:00:00"
RES_utr_extension="--cpus-per-task=2 --mem=64G --time=12:00:00"
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

# STARsolo: extra options, e.g. "--soloBarcodeReadLength 0" if GEX R1 is longer than 28 bp
STAR_EXTRA_OPTS="--soloBarcodeReadLength 0"
DEDUP_PARALLEL=4          # umi_tools chunks deduplicated in parallel (memory hungry)

# 3' UTR extension
UTR_EXT_BIN=100
UTR_EXT_DISTANCES=0,100,200,300,400

# batch-enriched cluster removal (batch_enriched_clusters_removal.R)
BATCH_FOLD=4                        # flag cluster if one batch > FOLD x another batch of the same stage group
LATE_STAGE_GROUPS="late"            # stage groups for which clusters are flagged when they contribute > LATE_MAX_FRACTION
LATE_MAX_FRACTION=0.8
MIN_WNN_FEATURES=200

# final cell type assignment
# optional table "cluster_number<TAB>lineage_label" for late clusters (see config/late_cluster_reassignment.tsv)
LATE_CLUSTER_TABLE=${PIPELINE_DIR}/config/late_cluster_reassignment.tsv
# cell types relabelled as 'unassigned' after inspection (paper: "Cxx"); space-separated, may be empty
UNASSIGN_CELL_TYPES="Cxx"

# IDR on pseudo-replicates (samples split with the 'idr_rep' column of samples.tsv)
IDR_FINAL_SET=0.05        # which IDR set builds the final peak set: 0.05 (paper) or 0.01
IDR_QC_HEATMAPS=yes

# optional motif analysis: MEME motif database (from https://meme-suite.org/meme/doc/download.html)
MEME_DB=${PROJECT_DIR}/motif_databases/CIS-BP_2.00/Caenorhabditis_elegans.meme
