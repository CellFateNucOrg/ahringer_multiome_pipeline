#!/usr/bin/env bash
#SBATCH --job-name=ja_multiome
#SBATCH --time=3-00:00:00                # walltime for the Nextflow *driver* (not the work)
#SBATCH --cpus-per-task=2
#SBATCH --mem=4G

# this script is not strictly necessary, as it just launches the job and exits immediately, but
# it does help to keep track of which stage was last run.
WORK_DIR=/mnt/meister.data/jsemple/20260728_ma_10x_multi_l3_PMW941_1/ahringer_pipeline_claude
stage=7
bash ${WORK_DIR}/run_pipeline.sh $stage


