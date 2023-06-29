#!/bin/bash
#
#SBATCH --output=/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses/data/fst/grenedalf_fst_%A.out.txt
#SBATCH --error=/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses/data/fst/grenedalf_fst_%A.err.txt
#SBATCH --time=23:59:00
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=40G

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${SLURM_JOBID}"

# run grenedalf
cd /Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses/data/fst/

# loop through each generation
for ii in {1..3}; do
echo -e "[$(date "+%Y-%m-%d %T")] generation = ${ii}"

grenedalf fst --window-type=genome --method=unbiased-nei \
--frequency-table-path "mergep_gen${ii}.csv" \
--pool-sizes "nflowers_gen${ii}.txt" \
--file-suffix "_gen${ii}"

done

echo -e "[$(date "+%Y-%m-%d %T")] JOB ID ${SLURM_JOBID} Done"
