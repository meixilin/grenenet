# Improved linear mixed model to test for correlations of allele frequency shifts with climates

## Math

To control for population structure, we used first 3 PCs from the PCA ($\bf{PC}$) calculated from the founder VCF. To reflect the changes in ecotype frequency over generations ($\Delta E$), the population structure matrix $\bf{C}$ is constructed as:

$$
\bf{C = \Delta E \times PC}
$$

By multiplying ecotype frequency changes with founder VCF PC, we decomposed two components that could impact the allele frequency changes ($\bf{\Delta P}$):

1.  selection signal
2.  population structure (reflected in ecotype frequency shifts)

Therefore, the equations are set up as following: 

$$
y = C\alpha + X\beta + Zu + \epsilon
$$

```
deltap ~ PC1 + PC2 + PC3 + bio1 + (1|site) 
```

## Scripts

1.  step0_prepare_data: reads in multiple data sources, and generate one `rdata` for downstream analyses
2.  step1_run_lmm: paraellelized worker script.
3.  step2_site_tau: calculate per site tau

## Submitting record

> 2023-09-09 15:20:45 2023-09-12 updates 2023-09-21 update the delta_ecotype tables

``` bash
WORKDIR='/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses'
RSCRIPT=${WORKDIR}/scripts/lmm_pc/step0_prepare_data.R
LOGFILE=${WORKDIR}/data/lmm_pc/logs/step0_prepare_data.log 
sbatchrr $RSCRIPT $LOGFILE # 418615
```

> 2023-09-10 17:16:00

``` bash
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

``` bash
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


# ENVS=("bio1" "bio2" "bio3" "bio4" "bio5" "bio6" "bio7" "bio8" "bio9" "bio10" "bio11" "bio12" "bio13" "bio14" "bio15" "bio16" "bio17" "bio18" "bio19" "tmin01" "tmin02" "tmin03" "tmin04" "tmin05" "tmin06" "tmin07" "tmin08" "tmin09" "tmin10" "tmin11" "tmin12" "tmax01" "tmax02" "tmax03" "tmax04" "tmax05" "tmax06" "tmax07" "tmax08" "tmax09" "tmax10" "tmax11" "tmax12" "tavg01" "tavg02" "tavg03" "tavg04" "tavg05" "tavg06" "tavg07" "tavg08" "tavg09" "tavg10" "tavg11" "tavg12" "prec01" "prec02" "prec03" "prec04" "prec05" "prec06" "prec07" "prec08" "prec09" "prec10" "prec11" "prec12" "srad01" "srad02" "srad03" "srad04" "srad05" "srad06" "srad07" "srad08" "srad09" "srad10" "srad11" "srad12" "wind01" "wind02" "wind03" "wind04" "wind05" "wind06" "wind07" "wind08" "wind09" "wind10" "wind11" "wind12" "vapr01" "vapr02" "vapr03" "vapr04" "vapr05" "vapr06" "vapr07" "vapr08" "vapr09" "vapr10" "vapr11" "vapr12" "longitude" "latitude" "altitude")

ENVS=("bio2" "bio3" "bio4" "bio5" "bio6" "bio7" "bio8" "bio9" "bio10" "bio13" "bio14" "bio15" "bio16" "bio17" "bio18" "tmin01" "tmin02" "tmin03" "tmin04" "tmin05" "tmin06" "tmin07" "tmin08" "tmin09" "tmin10" "tmin11" "tmin12" "tmax01" "tmax02" "tmax03" "tmax04" "tmax05" "tmax06" "tmax07" "tmax08" "tmax09" "tmax10" "tmax11" "tmax12" "tavg01" "tavg02" "tavg03" "tavg04" "tavg05" "tavg06" "tavg07" "tavg08" "tavg09" "tavg10" "tavg11" "tavg12" "prec01" "prec02" "prec03" "prec04" "prec05" "prec06" "prec07" "prec08" "prec09" "prec10" "prec11" "prec12" "srad01" "srad02" "srad03" "srad04" "srad05" "srad06" "srad07" "srad08" "srad09" "srad10" "srad11" "srad12" "wind01" "wind02" "wind03" "wind04" "wind05" "wind06" "wind07" "wind08" "wind09" "wind10" "wind11" "wind12" "vapr01" "vapr02" "vapr03" "vapr04" "vapr05" "vapr06" "vapr07" "vapr08" "vapr09" "vapr10" "vapr11" "vapr12" "longitude" "latitude" "altitude")
GENS=("1" "2" "3")
for envvar in "${ENVS[@]}"; do
for gen in "${GENS[@]}"; do
sleep 2
submit_lmm $envvar $gen
done
done
```
