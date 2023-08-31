# Local adaptation theory in ecotype frequency shifts

## Math

### Per ecotype

For a locally-adapted population/ecotype $i$, we can write its fitness at another location as the following.

$$
w = w_\text{max}\exp(-\frac{(x-x_0)^2}{V_s})
$$

We also know that for an ecotype, its ecotype frequency in the next generation $p_1$, compared with the starting ecotype frequency $p_0$, follows this relationship:

$$
p_1 = p_0\frac{w}{\bar{w}}
$$

Where $\bar{w}$ is the average fitness of all the ecotypes in experimental sites. We can divide $\bar{w}$ on both sides of the equation:

$$
\frac{w}{\bar{w}} = \frac{w_\text{max}}{\bar{w}}\exp(-\frac{(x-x_0)^2}{V_s})
$$

Given that we have $\frac{p_1}{p_0} = \frac{w}{\bar{w}}$,

$$
\frac{p_1}{p_0} = \frac{w_\text{max}}{\bar{w}}\exp(-\frac{(x-x_0)^2}{V_s})
$$

Taking the log

$$
\log(\frac{p_1}{p_0}) = \log(\frac{w_\text{max}}{\bar{w}}) - \frac{1}{V_s}(x-x_0)^2
$$

We can use linear regression to estimate $V_s$

$$
y = \log(\frac{p_1}{p_0}) \\
a = \log(\frac{w_\text{max}}{\bar{w}})  \\
b = -\frac{1}{V_s} \\
x = (x-x_0)^2
$$

### Per experimental site

When taking this theory to a broader level, and considering all 231 ecotypes as one population, across the 231 ecotypes, we now have the same $w_\text{max}$. Here, we are evaluating if the strength of selection at different experimental sites is variable. To avoid confusion, we write this $V_s$ equivalent as $\tau$.

$$
w = w_\text{max}\exp(-\frac{(x-x_0)^2}{\tau})
$$

## Scripts

1.  step0_prepare_data: the one reads in multiple data sources, and generate one `rdata` for downstream analyses
2.  step1_ecotype_Vs: calculate per ecotype Vs
3.  step2_site_tau: calculate per site tau

## Submitting record

> Wed Aug 30 12:46:57 2023


```bash
WORKDIR='/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses'
RSCRIPT=${WORKDIR}/scripts/local_adaptation/step1_ecotype_Vs.R
LOGFILE=${WORKDIR}/data/local_adaptation/step1_ecotype_Vs.log # 385372

RSCRIPT=${WORKDIR}/scripts/local_adaptation/step2_site_tau.R
LOGFILE=${WORKDIR}/data/local_adaptation/step2_site_tau.log
```



