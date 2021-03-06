\name{rrisk.samplesBgr}
\alias{rrisk.samplesBgr}
\title{Plot the Gelman-Rubin convergence statistic}
\usage{
  rrisk.samplesBgr(node, beg=samplesGetBeg(),
    end=samplesGetEnd(), firstChain=samplesGetFirstChain(),
    lastChain=samplesGetLastChain(), thin=samplesGetThin(),
    bins=50, plot=TRUE, ...)
}
\arguments{
  \item{node}{character vector of length 1, name of a
  variable in the model}

  \item{beg}{argument to select a slice of monitored values
  corresponding to iterations \code{beg:end}}

  \item{end}{argument to select a slice of monitored values
  corresponding to iterations \code{beg:end}}

  \item{firstChain}{argument to select a sub group of
  chains to calculate the Gelman-Rubin convergence
  statistics for. Number of chains must be larger than one}

  \item{lastChain}{argument to select a sub group of chains
  to calculate the Gelman-Rubin convergence statistics for.
  Number of chains must be larger than one}

  \item{thin}{only use every thin-th value of the stored
  sample for statistics}

  \item{bins}{number of blocks}

  \item{plot}{logical, whether to plot the BGR statistics
  or only return the values. If \code{TRUE}, values are
  returned invisibly}

  \item{...}{further graphical parameters as in \code{par}}
}
\description{
  This function calculates and plots the Gelman-Rubin
  convergence statistic, as modified by Brooks and Gelman
  (1998).
}
\details{
  This function is an alias of the function
  \code{\link{samplesBgr}} from the package \pkg{BRugs}.
  The original function was modified to create a diagnostic
  plots in other format as it was implemented in the
  original version. For more details see the function
  \code{\link{samplesBgr}} from the package \pkg{BRugs}.
  \cr \cr This function is not intended to be called
  directly but is internally called by
  \code{\link{diagnostics}} within the functions
  \code{\link{rrisk.BayesZIP}} and
  \code{\link{rrisk.BayesPEM}}.
}
\examples{
\dontrun{rrisk.samplesBgr("se")}
}
\keyword{manip}

