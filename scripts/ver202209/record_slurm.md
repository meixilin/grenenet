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

double check on the NaN scaled deltaP

```bash
# Wed Aug 31 16:26:01 2022
run_rscript step1.x_checkDeltaAF.R
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

# loop through all the vars
# JOBs: 101064 ~ 101122
# JOBID 101119 101118 101117 101116 101115 101114 101113 101112 101111 101110 101109 101107 101106 101105 101104 101103 101086 101085 101084 101083 101082 101081 101079 101078 101077 101076 101075 101074 101059 101122 101121 101120 101102 101101 101100 101099 101098 101097 101096 101095 101094 101093 101092 101090 101089 101088 101087
ENVARRAY=("LATITUDE" "LONGITUDE" "bio2"  "bio3"  "bio4"  "bio5"  "bio6"  "bio7"  "bio8"  "bio9"  "bio10" "bio11" "bio12" "bio13" "bio14" "bio15" "bio16" "bio17" "bio18" "bio19" "prec1"  "prec2"  "prec3"  "prec4"  "prec5"  "prec6"  "prec7"  "prec8"  "prec9"  "prec10" "prec11" "prec12" "tmin1"  "tmin2"  "tmin3"  "tmin4"  "tmin5"  "tmin6"  "tmin7"  "tmin8"  "tmin9"  "tmin10" "tmin11" "tmin12" "tmax1"  "tmax2"  "tmax3"  "tmax4"  "tmax5"  "tmax6"  "tmax7"  "tmax8"  "tmax9"  "tmax10" "tmax11" "tmax12")
for ii in "${ENVARRAY[@]}"; do
    sleep 5
    echo $ii
    sbatch_rscript ${MYRSCRIPT} "$ii"
done
```

## improve the comparisons

```bash
# Wed Aug 31 16:41:35 PDT 2022
sbatch_rscript() {
    local RSCRIPT=${1}
    local ENVVAR=${2}
    cd /Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/scripts/ver202209/
    sbatch --job-name=$(basename ${RSCRIPT/.R}) slurm_wrapper.sh ${RSCRIPT} ${ENVVAR}
}

MYRSCRIPT="/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/scripts/ver202209/fig3v2_LMregressionDeltaAF.R"
sbatch_rscript ${MYRSCRIPT} "bio1" # 101154
```



