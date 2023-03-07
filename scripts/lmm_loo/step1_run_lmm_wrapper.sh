#!/bin/bash
#
#SBATCH --output=/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses/logs/lmm_loo/step1_run_lmm/step1_run_lmm_%A_%a.out.txt
#SBATCH --error=/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses/logs/lmm_loo/step1_run_lmm/step1_run_lmm_%A_%a.err.txt
#SBATCH --time=23:59:00
#SBATCH --array=1-299
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=4G

sleep $((RANDOM % 20))
set -euo pipefail

HOMEDIR="/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses/"
COMMITID=$(git --git-dir="${HOMEDIR}/.git" --work-tree="${HOMEDIR}/" rev-parse master)
RSCRIPT=${HOMEDIR}/scripts/lmm_loo/step1_run_lmm.R

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${SLURM_JOBID} (${SLURM_ARRAY_JOB_ID}-${SLURM_ARRAY_TASK_ID}). Commit ID ${COMMITID} ..."

Rscript --vanilla ${RSCRIPT} ${SLURM_ARRAY_TASK_ID}

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${SLURM_JOBID} Done"

