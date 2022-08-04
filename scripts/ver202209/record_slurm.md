# step0.3: calculate distance 

```bash
sbatch_rscript() {
    local RSCRIPT=${1}
    sbatch --job-name=$(basename ${RSCRIPT/.R}) slurm_wrapper.sh ${RSCRIPT}
}

RSCRIPT="/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/scripts/ver202209/step0.3_dist_AF.R"
sbatch_rscript ${RSCRIPT}

```

