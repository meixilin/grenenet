#!/bin/bash
#SBATCH --output=/central/groups/carnegie_poc/meixilin/grenenet/analyses/logs/wrapper_averageAF.out.txt
#SBATCH --error=/central/groups/carnegie_poc/meixilin/grenenet/analyses/logs/wrapper_averageAF.err.txt
#SBATCH --time=23:59:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G

# @version         v0
# @script          sbatch wrapper_averageAF.sh
# @description     Wrapper to calculate average AF for all sites

# Author: Meixi Lin
# 2023-01-17 15:39:11

################################################################################
## import packages

set -eo pipefail

eval "$(/central/groups/carnegie_poc/meixilin/software/miniconda3/bin/conda shell.bash hook)"
conda activate grene

################################################################################
## def functions

################################################################################
## def variables

# bash variables
SITE=${1}
YEAR=${2}
CHROM=${3}
# directories
HOMEDIR=/central/groups/carnegie_poc/meixilin/grenenet/analyses
COMMITID=$(git --git-dir="${HOMEDIR}/.git" --work-tree="${HOMEDIR}" rev-parse master)
WORKDIR=${HOMEDIR}/data/averageAF

mkdir -p ${WORKDIR}
mkdir -p ${WORKDIR}/stats

# scripts
WORKSCRIPT=${HOMEDIR}/scripts/averageAF/grene_averageAF.py

################################################################################
## main
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${SLURM_JOBID}; git commit id: ${COMMITID}"

# run the scripts
python ${WORKSCRIPT} ${SITE} ${YEAR} ${CHROM}

conda deactivate
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${SLURM_JOBID} Done"

