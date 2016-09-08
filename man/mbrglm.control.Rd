\name{mbrglm.control}
\alias{mbrglm.control}
\title{
Auxiliary for Controlling MBRGLM Fitting
}
\description{
Auxiliary function as user interface for \code{\link{mbrglm}}
fitting. Typically only used when calling \code{\link{mbrglm}} or \code{\link{mbrglm.fit}}.
}
\usage{
mbrglm.control(mbr.epsilon = 1e-06, mbr.maxit = 500, 
  mbr.trace = FALSE, ...)
}

\arguments{
  \item{mbr.epsilon}{positive convergence tolerance for the iteration
    described in \code{\link{mbrglm.fit}}.}
  \item{mbr.maxit}{integer giving the maximum number of iterations for
  the iteration in \code{\link{mbrglm.fit}}.}
  \item{mbr.trace}{logical indicating if output should be produced for
    each iteration.}
  \item{\dots}{further arguments passed to or from other methods.}
}
\details{
The control argument of \code{\link{mbrglm}}is by default passed to the control argument of \code{\link{mbrglm.fit}}, which uses its elements as arguments to glm.control: the latter provides defaults.

If \code{mbr.trace=TRUE}, calls to \code{\link{cat}} produce the iteration number, the parameter estimates and the current value of the modified scores for each iteration.
}
\value{
A list with the arguments as components.
}
\references{
Kenne Pagui, E. C., Salvan, A. and Sartori, N. (2016). Median biasreduction of maximum likelihood estimates. \emph{http://arxiv.org/abs/1604.04768}.
}
\seealso{
  \code{\link{mbrglm.fit}}, the fitting procedure used by
  \code{\link{mbrglm}}.
}
\keyword{iteration}




