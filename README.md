# Welcome to the dev-version of the 

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
library(flip)
```

_A univariate analysis_

Testing the symmetry around 0 in a one sample (i.e. equivalent to one sample t-test) 

```r
set.seed(1)
y=rnorm(10)+.5
res=flip(y)
summary(res)
```

```
##  Call:
##  flip(Y = y) 
## 1023 permutations.
##   Test  Stat tail p-value sig.
## Y    t 2.561   ><  0.0293    *
```

* * *

## References

*Kenne Pagui, E. C., Salvan, A. and Sartori, N. (2016). Median bias reduction of maximum likelihood estimates. http://arxiv.org/abs/1604.04768.*
