\name{mbrglm}
\alias{mbrglm}
\alias{mbrglm.fit}
\alias{print.mbrglm}
\alias{summary.mbrglm}
\alias{print.summary.mbrglm}
\title{
Median bias reduction in Binomial-response GLMs
}
\description{
Fits binomial-response GLMs using the median bias-reduction method proposed in Kenne Pagui et al. (2016, Section 3). 
The proposed method is obtained by modifying the score equation in such a way that the solution is an approximately median 
unbiased estimator for each parameter component. The median bias-reduction method enjoys several good properties with respect to the maximum likelihood.
In particular,  the resulting estimator is component-wise median unbiased with and error of order (\eqn{\mathop{\rm
O}(n^{-1})}{O(n^{-1})}) and is equivariant under joint reparameterizations that transform each parameter component separately. 
It has the same asymptotic distribution as the maximum likelihood estimator. Moreover, the resulting estimates and their corresponding 
standard errors are always finite while the maximum likelihood estimates can be infinite in situations where complete or quasi separation occurs.
}
\usage{
mbrglm(formula, family = binomial, data, weights, subset, na.action, start = NULL, 
 etastart, mustart, offset, model = TRUE, method = "mbrglm.fit", x = FALSE, 
  y = TRUE, contrasts = NULL, control.glm = glm.control(), 
  control.mbrglm = mbrglm.control(), ...)

mbrglm.fit(x, y, weights = rep(1, nobs), start = NULL,  etastart = NULL, mustart = NULL,
 offset = rep(0, nobs), family = binomial(), control = glm.control(), 
 control.mbrglm = mbrglm.control(),  intercept = TRUE)
}
\arguments{
\item{formula}{an object of class \code{\link{formula}} (or one that can be coerced to that class): a symbolic description of the model to be fitted.}

\item{family}{a description of the error distribution and link function to be used in the model. For glm this can be a character string naming a family function, a family function or the result of a call to a family function. For mbrglm.fit only the third option is supported. (See \code{\link{family}} for details of family functions.) mbrglm currently supports only the "binomial" family with links "logit", "probit", "cloglog".}

\item{data}{an optional data frame, list or environment (or object coercible by \code{\link{as.data.frame}} to a data frame) containing the variables in the model. If not found in data, the variables are taken from environment(formula), typically the environment from which glm is called.}

\item{weights}{an optional vector of 'prior weights' to be used in the fitting process. Should be NULL or a numeric vector.}

\item{subset}{an optional vector specifying a subset of observations to be used in the fitting process.}

\item{na.action}{a function which indicates what should happen when the data contain NAs. The default is set by the na.action setting of \code{\link{options}}, and is \code{\link{na.fail}} if that is unset. The 'factory-fresh' default is \code{\link{na.omit}}. Another possible value is NULL, no action. Value \code{\link{na.exclude}}can be useful.}

\item{start}{starting values for the parameters in the linear predictor.}
\item{etastart}{starting values for the linear predictor.}

\item{mustart}{starting values for the vector of means.}

\item{offset}{this can be used to specify an a priori known component to be included in the linear predictor during fitting. This should be NULL or a numeric vector of length equal to the number of cases. One or more \code{\link{offset}} terms can be included in the formula instead or as well, and if more than one is specified their sum is used. See \code{\link{model.offset}}.}

\item{control}{a list of parameters for controlling the fitting process. For glm.fit this is passed to \code{\link{glm.control}}.}

\item{intercept}{logical. Should an intercept be included in the null model?}

\item{model}{a logical value indicating whether model frame should be included as a component of the returned value.}

\item{method}{the method to be used for fitting the model. The unique method is "mbrglm.fit", which uses the median modified score function to estimate the parameters.}

\item{x}{For mbrglm: logical values indicating whether the model matrix used in the fitting process should be returned as components of the returned value.}

\item{y}{For mbrglm: logical values indicating whether the response vector  used in the fitting process should be returned as components of the returned value.}

\item{contrasts}{an optional list. See the contrasts.arg of model.matrix.default.}

\item{control.glm}{\code{control.glm} replaces the \code{control}
  argument in \code{\link{glm}} but essentially does the same job. It
  is a list  of parameters to control \code{\link{glm.fit}}. See the
  documentation of  \code{glm.control1} for details.}

\item{control.mbrglm}{a list of parameters for controlling the fitting process when method="mbrglm.fit". See documentation \code{\link{mbrglm.control}} for details.}

\item{...}{additional arguments passed to or from other methods.}
}
\details{
\code{mbrglm.fit} is the workhorse function for fitting the model using
the median bias-reduction method.

The main iteration of \code{mbrglm.fit} consists to calculate the required quantities for the construction
of the modified iterative re-weighted least square which involves the modification term of the score function in the
adjusted dependent variable.

Iteration is repeated until either the iteration limit has been reached or the Euclidean distance  of the median modified scores is less than some specified positive constant (see the \code{mbr.maxit} and
\code{mbr.epsilon} arguments in \code{\link{mbrglm.control}}).
}
\value{
\code{\link{mbrglm}} returns an object of class \code{"mbrglm"}. A
  \code{"mbrglm"} object inherits first from \code{"glm"} and then from
  \code{"lm"} and is a list containing the following components:
  \item{coefficients}{a named vector of coefficients.}
  \item{residuals}{ Pearson's residual in the final iteration of the IWLS fit. Since cases with zero weights are omitted, their working residuals are NA.}
  \item{fitted.values}{the fitted mean values, obtained by transforming the linear predictors by the inverse of the link function.}
  \item{rank}{the numeric rank of the fitted linear model.}
  \item{family}{the \code{\link{family}} object used.}
  \item{linear.predictors}{the linear fit on link scale.}
  \item{deviance}{up to a constant, minus twice the maximized log-likelihood. Where sensible, the constant is chosen so that a saturated model has deviance zero.}
  \item{null.deviance}{The deviance for the null model, comparable with deviance. The null model will include the offset, and an intercept if there is one in the model. Note that this will be incorrect if the link function depends on the data other than through the fitted mean: specify a zero offset to force a correct calculation.}
  \item{weights}{the working weights, that is the weights in the final iteration of the IWLS fit.}
  \item{prior.weights}{the weights initially supplied, a vector of 1s if none were.}
  \item{df.residual}{the residual degrees of freedom.}
  \item{df.null}{the residual degrees of freedom for the null model.}
  \item{y}{if requested (the default) the y vector used. (It is a vector even for a binomial model.)}
  \item{x}{if requested, the model matrix.}
  \item{converged}{logical. Was the modified IWLS algorithm judged to have converged?}
  \item{boundary}{logical. Is the fitted value on the boundary of the attainable values?}
  \item{ModifiedScores}{the vector of the median modified scores for the
    parameters at the final iteration.}
  \item{FisherInfo}{the Fisher information matrix evaluated at the
    resulting estimates. Only available when \code{method = "mbrglm.fit"}.}
    \item{FisherInfoInvs}{the inverse of Fisher information matrix evaluated at the
    resulting estimates.}
    \item{nIter}{the number of iterations that were required until
    convergence. Only available when \code{method = "mbrglm.fit"}.}
    \item{model}{if requested (the default), the model frame.}
    \item{formula}{the formula supplied.}
    \item{terms}{the terms object used.}
    \item{data}{the data argument.}
    \item{offset}{the offset vector used.}
    \item{control.mbrglm}{the \code{control.mbrglm} argument that was passed to
  \code{mbrglm}. Only available when \code{method = "mbrglm.fit"}.}
   \item{contrasts}{(where relevant) the contrasts used.}
}
\note{
  1.  'mbrglm' and 'mbrglm.fit' were written using as basis structure
 the code of  'brglm' and 'brglm.fit', respectively. The functions 'brglm' and 'brglm.fit' are implemented in the R package brglm version 0.5-9 by Ioannis Kosmidis. While,  'print.mbrglm', 'summary.mbrglm' and  'print.summary.mbrglm' are  modifications of  'print.glm', 'summary.glm' and  'print.summary.glm', respectively.
     }

\references{
Kenne Pagui, E. C., Salvan, A. and Sartori, N. (2016). Median bias
reduction of maximum likelihood estimates. \emph{http://arxiv.org/abs/1604.04768}.
}

\seealso{brglm, brglm.fit, \code{\link{glm}}, \code{\link{glm.fit}}}

\examples{
 ## First example
library(brglm)     
data(endo)
# Fit the GLM using maximum likelihood
endo.glm <- glm(HG~NV+PI+EH,family=binomial,data=endo)
## Mean bias-reduced fit
endo.brglm<-brglm(HG~NV+PI+EH,family=binomial,data=endo)
## Median bias-reduced fit
endo.mbrglm<-mbrglm(HG~NV+PI+EH,family=binomial,data=endo)
endo.glm
endo.brglm
endo.mbrglm

# Now other links
update(endo.mbrglm, family = binomial(probit))
update(endo.mbrglm, family = binomial(cloglog))

##------------------------
## paper by Andrey Gelman et al. 2008. Annals of applied Statistics.
## application to binomial
## example 4.2
##----------------------

# first way

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

# in alternative
fit.glm<-glm(cbind(z,weights-z)~z.x,family=binomial)
fit.brglm<-brglm(cbind(z,weights-z)~z.x,family=binomial)
fit.mbrglm<-mbrglm(cbind(z,weights-z)~z.x,family=binomial)
fit.glm
fit.brglm
fit.mbrglm

##----------------------------------------
# Rasch model: 100 subjects and 5 items
##----------------------------------------

I <- 5
S <-  100

## function to generate data
gendata.M <- function(gamma, alpha, beta)
{
  I <- length(alpha) 
  S <- length(gamma) 
  data.y <- matrix(0, nrow=S, ncol=I)
  for(i in 1:I)
  {
    mui <- plogis(alpha[i] + gamma * beta[i])
    data.y[,i] <- rbinom(S, size=1, prob=mui)
  }
  return(data.y)      
}

alphas <- c(0.0,  0.7,  1.6,  0.6, -0.5)
betas <- rep(1,I)
gammas <- rnorm(S)

y <- gendata.M(gammas,alphas,betas)

y.dat <- data.frame(y=y[1:(S*I)],subject=factor(rep(1:S,I)),item=factor(rep(1:I,each=S)))

\dontrun{
fit.glm <- glm(y~subject-1+item,family=binomial,data=y.dat)
fit.brglm <- brglm(y~subject-1+item,family=binomial,data=y.dat)
fit.mbrglm <- mbrglm(y~subject-1+item,family=binomial,data=y.dat)
}
summary(fit.glm)
summary(fit.brglm)
summary(fit.mbrglm)

}
\keyword{models}
\keyword{regression}
\keyword{iteration}