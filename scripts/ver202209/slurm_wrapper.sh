#!/bin/bash
#
##SBATCH --job-name=xx
#SBATCH --output=/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/logs/ver202209/slurm_wrapper.out.txt
#SBATCH --error=/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/logs/ver202209/slurm_wrapper.err.txt
#SBATCH --time=23:59:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=20G

# to submit:
# sbatch_rscript() {
#     local RSCRIPT=${1}
#     sbatch --job-name=$(basename ${RSCRIPT/.R}) slurm_wrapper.sh ${RSCRIPT}
# }

set -eo pipefail

RSCRIPT=${1}
INPUT=${2}
RSCRIPTNAME=$(basename ${RSCRIPT/.R})

if [ ! -f ${RSCRIPT} ]; then
    echo -e "[$(date "+%Y-%m-%d %T")] FAIL File not found ${RSCRIPT}" 1>&2
    # exit 1
fi

HOMEDIR="/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/"
COMMITID=$(git --git-dir="${HOMEDIR}/.git" --work-tree="${HOMEDIR}/" rev-parse master)

cd ${HOMEDIR}

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${SLURM_JOBID}; GIT commit id ${COMMITID}; Running ${RSCRIPT} ${INPUT}..."

Rscript --vanilla ${RSCRIPT} ${INPUT} &> "logs/ver202209/${RSCRIPTNAME}.log"

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${SLURM_JOBID} Done"

