# step0.3: calculate distance 

```bash
sbatch_rscript() {
    local RSCRIPT=${1}
    cd /Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/scripts/ver202209/
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

run_rscript "fig2.2_nmdsSURF_af.R"
```

# fig3 AF by environmental data

## first load gff files and find the gene locations

```bash
run_rscript "step0.4_loadGFF.R"
```


## load the delta AF files

```bash
sbatch_rscript "/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/scripts/ver202209/step1.0_loadDeltaAF.R" # 100586
```

## run the comparisons

> Wed Aug 31 01:22:31 PDT 2022

```bash
sbatch_rscript() {
    local RSCRIPT=${1}
    local ENVVAR=${2}
    cd /Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/scripts/ver202209/
    sbatch --job-name=$(basename ${RSCRIPT/.R}) slurm_wrapper.sh ${RSCRIPT} ${ENVVAR}
}

MYRSCRIPT="/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/scripts/ver202209/fig3_LMregressionDeltaAF.R"
sbatch_rscript ${MYRSCRIPT} "bio1" # 100851

```





