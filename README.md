# Welcome to the ongoing work of the library mbrglm

# [library mbrglm]

* * *

## Set up

To **install** this github version type (in R):

    #if devtools is not installed yet: 
    # install.packages("devtools") 
    library(devtools)
    install_github("livioivil/mbrglm")


* * *

## Some examples

```r
library(mbrglm)
data(endo)
library(brglm)
```

```r
# Fit the GLM using maximum likelihood
endo.glm <- glm(HG~NV+PI+EH,family=binomial,data=endo)
## Mean bias-reduced fit usind the R package brglm
endo.brglm<-brglm(HG~NV+PI+EH,family=binomial,data=endo)
## Median bias-reduced fit
endo.mbrglm<-mbrglm(HG~NV+PI+EH,family=binomial,data=endo)
endo.glm
endo.brglm
endo.mbrglm
```

```
endo.glm

## 
## Call:  glm(formula = HG ~ NV + PI + EH, family = binomial, data = endo)
## 
## Coefficients:
## (Intercept)           NV           PI           EH  
##     4.30452     18.18556     -0.04218     -2.90261  
## 
## Degrees of Freedom: 78 Total (i.e. Null);  75 Residual
## Null Deviance:       104.9 
## Residual Deviance: 55.39     AIC: 63.39

endo.brglm

## 
## Call:  brglm(formula = HG ~ NV + PI + EH, family = binomial, data = endo) 
## 
## Coefficients:
## (Intercept)           NV           PI           EH  
##     3.77456      2.92927     -0.03475     -2.60416  
## 
## Degrees of Freedom: 78 Total (i.e. Null);  75 Residual
## Deviance:        56.5754 
## Penalized Deviance: 48.0745  AIC: 64.5754

endo.mbrglm

## 
## Call:  mbrglm(formula = HG ~ NV + PI + EH, family = binomial, data = endo) 
## 
## Standardized Pearson residual:
##     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.  
## -1.52800  -0.49990  -0.23260   0.05676   0.15680   5.80200  
## 
## Coefficients:
## (Intercept)           NV           PI           EH  
##     3.96936      3.86921     -0.03868     -2.70793  
## 
## Degrees of Freedom: 78 Total (i.e. Null);  75 Residual
## Deviance:        55.8679
```


* * *

## References

*Kenne Pagui, E. C., Salvan, A. and Sartori, N. (2016). Median bias reduction of maximum likelihood estimates. http://arxiv.org/abs/1604.04768.*
