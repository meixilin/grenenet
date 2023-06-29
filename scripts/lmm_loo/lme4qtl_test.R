# devtools::install_github("variani/lme4qtl")
library(lme4qtl)

data(dat40, package = "lme4qtl")

dim(dat40)
#> [1] 234  10
dim(kin2)
#> [1] 234 234

head(dat40)
#>     ID  trait1  trait2 AGE FAMID  FA  MO SEX trait1bin trait2bin
#> 7  101 8.41954 9.67925  50    10   0   0   1         0         0
#> 14 102 5.47141 4.31886  44    10   0   0   2         0         0
#> 21 103 9.66547 7.00735  34    10 101 102   2         0         0
#> 28 104 6.27092 8.59257  41    10 101 102   1         0         0
#> 35 105 7.96814 7.60801  36    10 101 102   1         0         0
#> 42 106 8.29865 8.17634  37    10 101 102   2         0         0
kin2[1:5, 1:5] # nuclear family with 2 parents and 3 kids
#> 5 x 5 sparse Matrix of class "dsCMatrix"
#>     11  12  13  14  15
#> 11 1.0 .   0.5 0.5 0.5
#> 12 .   1.0 0.5 0.5 0.5
#> 13 0.5 0.5 1.0 0.5 0.5
#> 14 0.5 0.5 0.5 1.0 0.5
#> 15 0.5 0.5 0.5 0.5 1.0
#> 

# does it matter for kin2 to be dsCMatrix
kin3 = as.matrix(kin2)

m1 <- relmatLmer(trait1 ~ AGE + SEX + (1|FAMID) + (1|ID), dat40, relmat = list(ID = kin2))
m1 <- relmatLmer(trait1 ~ AGE + SEX + (1|FAMID) + (1|ID), dat40, relmat = list(ID = kin3))
#> boundary (singular) fit: see ?isSingular
#> boundary (singular) fit: see ?isSingular
m1
#> Linear mixed model fit by REML ['lmerMod']
#> Formula: trait1 ~ AGE + SEX + (1 | FAMID) + (1 | ID)
#>    Data: dat40
#> REML criterion at convergence: 963.3853
#> Random effects:
#>  Groups   Name        Std.Dev.
#>  ID       (Intercept) 2.2988  
#>  FAMID    (Intercept) 0.0000  
#>  Residual             0.7856  
#> Number of obs: 224, groups:  ID, 224; FAMID, 39
#> Fixed Effects:
#> (Intercept)          AGE         SEX2  
#>    7.563248     0.008314    -0.364197  
#> convergence code 0; 1 optimizer warnings; 0 lme4 warnings
#> 

library(lme4)

# try run the model using lme4
m2 <- lme4::lmer(trait1 ~ AGE + SEX + (1|FAMID) + (1|ID), dat40, 
                 control=lmerControl(check.nobs.vs.nlev = "ignore",
                                     check.nobs.vs.rankZ = "ignore",
                                     check.nobs.vs.nRE="ignore"))
m2
# Linear mixed model fit by REML ['lmerMod']
# Formula: trait1 ~ AGE + SEX + (1 | FAMID) + (1 | ID)
# Data: dat40
# REML criterion at convergence: 986.6654
# Random effects:
#     Groups   Name        Std.Dev.
# ID       (Intercept) 1.459   
# FAMID    (Intercept) 1.283   
# Residual             1.277   
# Number of obs: 224, groups:  ID, 224; FAMID, 39
# Fixed Effects:
#     (Intercept)          AGE         SEX2  
# 7.34300      0.01402     -0.49361 

# check if the model's residual matrix is as specified



