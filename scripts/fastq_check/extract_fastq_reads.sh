#!/bin/bash
#

# @version         v0
# @script          bash extract_fastq_reads.sh
# @description     extract the fastqc output on 'Total sequences'

# Author: Meixi Lin
# 2023-01-26 09:50:23

################################################################################
## import packages

set -eo pipefail

################################################################################
## def functions

################################################################################
## def variables
# directories
HOMEDIR=/central/groups/carnegie_poc/meixilin/grenenet/analyses
FASTQCDIR=/central/groups/carnegie_poc/lczech/grenephase1/qc/fastqc
COMMITID=$(git --git-dir="${HOMEDIR}/.git" --work-tree="${HOMEDIR}" rev-parse master)

WORKDIR=${HOMEDIR}/data/fastq_check
mkdir -p ${WORKDIR}
OUTFILE='fastqc_totalseq.txt'

################################################################################
## main
echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${SLURM_JOBID}; git commit id: ${COMMITID}"
cd $FASTQCDIR

for FILE in *.html; do
NREADS=$(tail -1 $FILE | sed 's/.*Total Sequences<\/td><td>//; s/<\/td>.*//')
echo -e "${FILE/_fastqc.html}\t${NREADS}" >> ${WORKDIR}/${OUTFILE}
done

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${SLURM_JOBID} Done"
