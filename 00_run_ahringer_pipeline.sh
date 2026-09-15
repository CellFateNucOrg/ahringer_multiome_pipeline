#!/usr/bin/env bash
#SBATCH --job-name=ja_multiome
#SBATCH --time=3-00:00:00                # walltime for the Nextflow *driver* (not the work)
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G


WORK_DIR=/mnt/meister.data/jsemple/20260728_ma_10x_multi_PMW941_1/ahringer_pipeline_claude
stage=2
${WORK_DIR}/run_pipeline.sh $stage


