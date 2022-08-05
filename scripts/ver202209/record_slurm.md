# step0.3: calculate distance 

```bash
sbatch_rscript() {
    local RSCRIPT=${1}
    sbatch --job-name=$(basename ${RSCRIPT/.R}) slurm_wrapper.sh ${RSCRIPT}
}

RSCRIPT="/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/scripts/ver202209/step0.3_dist_AF.R"
sbatch_rscript ${RSCRIPT} # 95478

```

# fig2.1 nmds

```bash
run_rscript() {
    local RSCRIPT=${1}
    local RSCRIPTNAME=$(basename ${RSCRIPT/.R})
    Rscript --vanilla ${RSCRIPT} &> "/home/mlin/safedata/meixilin/grenenet/logs/ver202209/${RSCRIPTNAME}.log"
    echo $?
}

run_rscript "fig2.1_nmds_af.R"
```

