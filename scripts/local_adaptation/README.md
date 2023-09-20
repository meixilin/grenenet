# Local adaptation theory in ecotype frequency shifts

## Math

### Per ecotype

For a locally-adapted population/ecotype $i$, we can write its fitness at another location $j$ as the following.

$$
w = w_\text{max}\exp(-\frac{(x-x_0)^2}{V_s})
$$

We also know that for an ecotype, its ecotype frequency in the next generation $p_1$, compared with the starting ecotype frequency $p_0$, follows this relationship:

$$
p_1 = p_0\frac{w}{\bar{w_j}}
$$

Where $\bar{w_j}$ is the average fitness of all the ecotypes in experimental site $j$. We can divide $\bar{w_j}$ on both sides of the equation:

$$
\frac{w}{\bar{w_j}} = \frac{w_\text{max}}{\bar{w_j}}\exp(-\frac{(x-x_0)^2}{V_s})
$$

Given that we have $\frac{p_1}{p_0} = \frac{w}{\bar{w_j}}$,

$$
\frac{p_1}{p_0} = \frac{w_\text{max}}{\bar{w_j}}\exp(-\frac{(x-x_0)^2}{V_s})
$$

Taking the log

$$
\log(\frac{p_1}{p_0}) = \log(\frac{w_\text{max}}{\bar{w_j}}) - \frac{1}{V_s}(x-x_0)^2
$$

We can use linear mixed model regression to estimate $V_s$, which are assumed to be constant across sites within ecotypes but variable between ecotypes.

$$
y = \log(\frac{p_1}{p_0}) \\
a' = \log(\frac{w_\text{max}}{\bar{w_j}})  \\
b = -\frac{1}{V_s} \\
x = (x-x_0)^2
$$

As $\bar{w_j}$ is variable across sites, but $V_s$ is assumed to be constant, we can set a random intercept model, where, $a'$ is now written in two parts:

$$
y = (u + a) + bx + \epsilon
$$

### Per experimental site

When taking this theory to a broader level, and considering all 231 ecotypes as one population, across the 231 ecotypes, we now have the same $w_\text{max}$. Here, we are evaluating if the strength of selection at different experimental sites is variable. To avoid confusion, we write this $V_s$ equivalent as $\tau$. $\tau$ is assumed to be constant across ecotypes within a site but variable between sites.

An advantage for this model is that the $\bar{w}$ is the same since the average fitness at the same site is the same across the ecotype evaluated.

$$
w = w_\text{max}\exp(-\frac{(x-x_0)^2}{\tau})
$$

For this model, we can use the regular linear regression model without the need for invoking the random intercept terms.

## Scripts

1.  step0_prepare_data: the one reads in multiple data sources, and generate one `rdata` for downstream analyses
2.  step1_ecotype_Vs: calculate per ecotype Vs
3.  step2_site_tau: calculate per site tau

## Submitting record

> ARCHIVE: 2023-08-30; 2023-09-13; 2023-09-15

```bash
# Mon Sep 18 18:31:16 PDT 2023
WORKDIR='/Carnegie/DPB/Data/Shared/Labs/Moi/Everyone/meixilin/grenenet/analyses'

RSCRIPT=${WORKDIR}/scripts/local_adaptation/step0_prepare_data.R
LOGFILE=${WORKDIR}/data/local_adaptation/logs/step0_prepare_data.log # 415354

RSCRIPT=${WORKDIR}/scripts/local_adaptation/step1_ecotype_Vs.R
LOGFILE=${WORKDIR}/data/local_adaptation/logs/step1_ecotype_Vs.log # 415356

RSCRIPT=${WORKDIR}/scripts/local_adaptation/step2_site_tau.R
LOGFILE=${WORKDIR}/data/local_adaptation/logs/step2_site_tau.log # 415357

RSCRIPT=${WORKDIR}/scripts/local_adaptation/plot_local_adaptation.R
LOGFILE=${WORKDIR}/data/local_adaptation/logs/plot_local_adaptation.log # 418601

sbatchrr ${RSCRIPT} ${LOGFILE}
```

