# step0.3: calculate distance 

```bash
sbatch_rscript() {
    local RSCRIPT=${1}
    cd /Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/scripts/ver202209/
    sbatch --job-name=$(basename ${RSCRIPT/.R}) slurm_wrapper.sh ${RSCRIPT} ${2} ${3}
}

RSCRIPT="/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/scripts/ver202209/step0.3_dist_AF.R"
sbatch_rscript ${RSCRIPT} # 95478

```

## load the delta AF files

```bash
sbatch_rscript "/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/scripts/ver202209/step1.0_loadDeltaAF.R" # 100586
```

## run the comparisons
100609 100616

```bash
sbatch_rscript \
"/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/scripts/ver202209/step1.1_regressionDeltaAF.R" \
"/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/data/AF/ver202209/haplotype/DeltaP/delta_freq_uniq.rds" \
"LATITUDE"

sbatch_rscript \
"/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/scripts/ver202209/step1.1_regressionDeltaAF.R" \
"/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/data/AF/ver202209/haplotype/DeltaP/scaled_delta_freq_uniq.rds" \
"bio1"

sbatch_rscript \
"/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/scripts/ver202209/step1.1_regressionDeltaAF.R" \
"/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/data/AF/ver202209/haplotype/DeltaP/scaled_delta_freq_uniq.rds" \
"bio4"

sbatch_rscript \
"/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/scripts/ver202209/step1.1_regressionDeltaAF.R" \
"/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/data/AF/ver202209/haplotype/DeltaP/scaled_delta_freq_uniq.rds" \
"bio12"

sbatch_rscript \
"/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/scripts/ver202209/step1.1_regressionDeltaAF.R" \
"/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/data/AF/ver202209/haplotype/DeltaP/scaled_delta_freq_uniq.rds" \
"LATITUDE"
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




