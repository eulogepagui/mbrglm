\name{endo}
\alias{endo}
\docType{data}
\title{
Endometrial cancer grade
}
\description{Study of 79 primarily diagnosed cases of endometrial cancer in order to explain the histology of the endometrium by putative risk factors. The \code{endo} data frame has 79 rows and 4 columns. Variable histology, \code{HG}. Three risk factors:  variable neovasculization, \code{NV},  variable pulsatility index, \code{PI} and variable endometrium height, \code{EH}.
}
\usage{data("endo")}
\format{
  A data frame with 79 observations on the following 4 variables.
  \describe{
    \item{\code{NV}}{a factor with two levels 1 (present) and 0 (absent) }
    \item{\code{PI}}{pulsatility index of arteria uterina}
    \item{\code{EH}}{endometrium height}
    \item{\code{HG}}{a factor with two levels classified as either 0 (grading 0-II) or 1 (grading III-IV)}
  }
}
\source{
Agresti, A. (2015).  \emph{Foundations of Linear and Generalized Linear Models}. John Wiley & Sons.
}
\examples{
data(endo)
glm(HG~NV+PI+EH,family=binomial,data=endo)
mbrglm(HG~NV+PI+EH,family=binomial,data=endo)
}
\keyword{datasets}
