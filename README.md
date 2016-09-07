# Welcome to the ongoing work of the library mbrglm

# [library mbrglm]

* * *

## Set up

To **install** this github version type (in R):

    #if devtools is not installed yet: 
    # install.packages("devtools") 
    library(devtools)
    install_github("eulogepagui/mbrglm")


* * *

## Some examples
**First example: endometrial dataset**

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
Now other links

```r
update(endo.mbrglm, family = binomial(probit))
update(endo.mbrglm, family = binomial(cloglog))
```
```
update(endo.mbrglm, family = binomial(probit))
## 
## Call:  mbrglm(formula = HG ~ NV + PI + EH, family = binomial(probit),  data = endo) 
## 
## Standardized Pearson residual:
##     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.  
## -1.47200  -0.52520  -0.26210   0.04931   0.17090   5.95500  
## 
## Coefficients:
## (Intercept)           NV           PI           EH  
##     1.98426      1.97083     -0.01661     -1.42457  
## 
## Degrees of Freedom: 78 Total (i.e. Null);  75 Residual
## Deviance:        57.0296
update(endo.mbrglm, family = binomial(cloglog))
## 
## Call:  mbrglm(formula = HG ~ NV + PI + EH, family = binomial(cloglog), data = endo) 
## 
## Standardized Pearson residual:
##     Min.   1st Qu.    Median      Mean   3rd Qu.      Max.  
## -1.62200  -0.48640  -0.24340   0.03919   0.08179   5.28400  
## 
## Coefficients:
## (Intercept)           NV           PI           EH  
##     3.11967      1.80369     -0.03714     -2.32510  
## 
## Degrees of Freedom: 78 Total (i.e. Null);  75 Residual
## Deviance:        54.8008
```
**second example: example 4.2 of paper by Andrey Gelman et al. 2008. Annals of applied Statistics.**

```r
x<-c(-0.86,-0.30,-0.05,0.73)
z.x<- (1/sqrt(4))*(x-mean(x))/sqrt(var(x))
weights<-rep(5,4)
z<-c(0,1,3,5)
y=z/weights
fit.glm<-glm(y~z.x,family=binomial,weights=weights)
fit.brglm<-brglm(y~z.x,family=binomial,weights=weights)
fit.mbrglm<-mbrglm(y~z.x,family=binomial,weights=weights)
fit.glm
fit.brglm
fit.mbrglm
```
```
fit.glm
## 
## Call:  glm(formula = y ~ z.x, family = binomial, weights = weights)
## 
## Coefficients:
## (Intercept)          z.x  
##    -0.08328     10.23079  
## 
## Degrees of Freedom: 3 Total (i.e. Null);  2 Residual
## Null Deviance:       15.79 
## Residual Deviance: 0.05474   AIC: 7.965

fit.brglm

## 
## Call:  brglm(formula = y ~ z.x, family = binomial, weights = weights) 
## 
## Coefficients:
## (Intercept)          z.x  
##     -0.2032       5.4207  
## 
## Degrees of Freedom: 3 Total (i.e. Null);  2 Residual
## Deviance:        1.0484 
## Penalized Deviance: 1.9875   AIC: 8.9585

fit.mbrglm

## 
## Call:  mbrglm(formula = y ~ z.x, family = binomial, weights = weights) 
## 
## Standardized Pearson residual:
##      Min.    1st Qu.     Median       Mean    3rd Qu.       Max.  
## -0.300500  -0.267100   0.007402   0.015900   0.290400   0.349300  
## 
## Coefficients:
## (Intercept)          z.x  
##     -0.1712       7.5715  
## 
## Degrees of Freedom: 3 Total (i.e. Null);  2 Residual
## Deviance:        0.2755
```
* * *

## References

*Kenne Pagui, E. C., Salvan, A. and Sartori, N. (2016). Median bias reduction of maximum likelihood estimates. http://arxiv.org/abs/1604.04768.*
