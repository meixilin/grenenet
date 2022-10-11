#!/bin/bash
#
##SBATCH --job-name=release025_check
#SBATCH --output=/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/logs/metadata_check/release025_check.out.txt
#SBATCH --error=/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/logs/metadata_check/release025_check.err.txt
#SBATCH --time=08:00:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=8G

set -eo pipefail

HOMEDIR="/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/"
COMMITID=$(git --git-dir="${HOMEDIR}/.git" --work-tree="${HOMEDIR}/" rev-parse master)
DUPLIST02="${HOMEDIR}/scripts/metadata_check/duplicated_reads_release02.csv"
DUPLIST05="${HOMEDIR}/scripts/metadata_check/duplicated_reads_release05.csv"
cd ${HOMEDIR}/scripts/metadata_check

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${SLURM_JOBID}; GIT commit id ${COMMITID}; Running ..."

while read p; do
  md5sum "$p" >> duplicated_reads_release02.md5sum.txt
done < ${DUPLIST02}

while read p; do
  md5sum "$p" >> duplicated_reads_release05.md5sum.txt
done < ${DUPLIST05}

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${SLURM_JOBID} Done"

