% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/modcomp.R
\name{modelcomp}
\alias{modelcomp}
\title{A Function for Running Sequential Mixed-Effects Model Comparisons}
\usage{
modelcomp(origmodel, resdf = data.frame(), compdf = data.frame(),
  bootstr = FALSE)
}
\arguments{
\item{bootstr}{should bootstrapping be used. Defaults to FALSE.}
}
\description{
This function allows you run mixed-effects model comparisons with bootstrapping switched on or off
}
\examples{
modelcomp()
}
\keyword{bootstrapping,}
\keyword{comparisons,}
\keyword{mixed}
