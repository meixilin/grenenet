# Improved linear mixed model to test for correlations of allele frequency shifts with climates

## Math

To control for population structure, we used first 10 PCs from the PCA ($\bf{PC}$) calculated from the founder VCF. To reflect the changes in ecotype frequency over generations ($\Delta E$), the population structure matrix $\bf{C}$ is constructed as:

$$
\bf{C = \Delta E \times PC}
$$

By multiplying ecotype frequency changes with founder VCF PC, we decomposed two components that could impact the allele frequency changes ($\Delta P$):

1.  selection signal
2.  population structure (reflected in ecotype frequency shifts)

## Scripts

1.  step0_prepare_data: reads in multiple data sources, and generate one `rdata` for downstream analyses
2.  step1_run_lmm: paraellelized worker script.
3.  step2_site_tau: calculate per site tau

## Submitting record

> 2023-09-09 15:20:45 2023-09-12 updates 2023-09-21 update the delta_ecotype tables

```bash
WORKDIR='/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses'
RSCRIPT=${WORKDIR}/scripts/lmm_pc/step0_prepare_data.R
LOGFILE=${WORKDIR}/data/lmm_pc/logs/step0_prepare_data.log 
sbatchrr $RSCRIPT $LOGFILE # 418615
```

> 2023-09-10 17:16:00

```bash
# for now test whether using PCs works and how many works
submit_lmm(){
    envvar=${1}
    mygen=${2}
    local WORKDIR='/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses'
    local RSCRIPT=${WORKDIR}/scripts/lmm_pc/step1_run_lmm.R
    local LOGFILE=${WORKDIR}/data/lmm_pc/logs/step1_run_lmm_${envvar}_gen${mygen}.log
    sbatchrr ${RSCRIPT} ${LOGFILE} ${envvar} ${mygen} 
}

# these were ran with code currently named as `ARCHIVE/step1_run_lmm.R`
submit_lmm 'bio1' '1'
submit_lmm 'bio1' '2'
submit_lmm 'bio1' '3'
submit_lmm 'bio11' '1'
submit_lmm 'bio11' '2'
submit_lmm 'bio11' '3'
submit_lmm 'bio19' '1'
submit_lmm 'bio19' '2'
submit_lmm 'bio19' '3'

# new version with precalculated and LD pruned deltap
submit_lmm 'bio19' '1' # 412313
```

> final version: using 5 PCs and LD pruned deltap

```bash
# Thu Sep 21 19:31:24 PDT 2023
submit_lmm 'bio1' '1' # 418628 
submit_lmm 'bio1' '2'
submit_lmm 'bio1' '3'
submit_lmm 'bio11' '1'
submit_lmm 'bio11' '2'
submit_lmm 'bio11' '3'
submit_lmm 'bio12' '1'
submit_lmm 'bio12' '2'
submit_lmm 'bio12' '3'
submit_lmm 'bio19' '1'
submit_lmm 'bio19' '2'
submit_lmm 'bio19' '3'
```


