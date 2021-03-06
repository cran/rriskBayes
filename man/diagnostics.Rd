\name{diagnostics}
\alias{diagnostics}
\title{Diagnostic plots for MCMC models provided by rrisk.BayesPEM and rrisk.BayesZIP functions}
\usage{
  diagnostics(nodes, plots=FALSE)
}
\arguments{
  \item{nodes}{character string, the name of parameters(s)}

  \item{plots}{logical, if \code{TRUE} the diagnostic plots
  will be displayed in separate windows}
}
\value{
  Returns \code{TRUE} if the user confirms convergence.
  Otherwise the function returns \code{FALSE}.
}
\description{
  This function provides a GUI for diagnostic plots to
  check convergence in Markov chain Monte-Carlo (MCMC)
  models provided by \code{\link{rrisk.BayesZIP}} and
  \code{\link{rrisk.BayesPEM}}.
}
\details{
  The argument \code{nodes} denotes the node(s) to be used
  for diagnostic plots of the MCMC chains. The user is
  interactively requested to confirm whether the
  convergence has been reached. In this case the function
  returns \code{TRUE} otherwise \code{FALSE}. \cr \cr This
  function is not intended to be called directly but is
  internally called within \code{\link{rrisk.BayesZIP}} or
  \code{\link{rrisk.BayesPEM}}.
}
\examples{
\dontrun{diagnostics(nodes)}
}
\seealso{
  For more details see documentation to the functions
  \code{\link{samplesBgr}}, \code{\link{plotHistory}} and
  \code{\link{plotDensity}} from the package \pkg{BRugs}.
}
\keyword{manip}

