\name{PEMGUI}
\alias{PEMGUI}
\title{GUI for Bayesian Prevalence estimation under misclassification (PEM)}
\usage{
  PEMGUI(x=20, n=20, k=10, prior.pi=c(1,19),
    prior.se=c(1,1), prior.sp=c(1,1), chains=3, burn=1000,
    update=10000, thin=1)
}
\arguments{
  \item{x}{scalar value for number of pools (\code{k>1}) or
  single outcomes (\code{k=1}) with positive test result}

  \item{n}{scalar value for number of pools tested
  (\code{k>1}) or the sample size in application study
  (\code{k=1})}

  \item{k}{scalar value for number of individual samples
  physically combined into one pool}

  \item{prior.pi}{numeric vector containing parameters of a
  beta distribution as prior for prevalence \code{pi}, e.g.
  \code{pi} \eqn{\sim} \code{prior.pi(*,*)=beta(*,*)}}

  \item{prior.se}{numeric vector containing parameters of a
  beta distribution as prior for sensitivity \code{se},
  e.g. \code{se} \eqn{\sim} \code{prior.se(*,*)=beta(*,*)}}

  \item{prior.sp}{numeric vector containing parameters of a
  beta distribution as prior for specificity \code{sp},
  e.g. \code{sp} \eqn{\sim} \code{prior.sp(*,*)=beta(*,*)}}

  \item{chains}{positive single numeric value, number of
  independent MCMC chains (default 3)}

  \item{burn}{positive single numeric value, length of the
  burn-in period (default 1000)}

  \item{update}{positive single numeric value, length of
  update iterations for estimation (default 10000)}

  \item{thin}{positive single numeric value (default 1).
  The samples from every kth iteration will be used for
  inference, where k is the value of thin. Setting
  \code{thin > 1} can help to reduce the autocorrelation in
  the sample.}
}
\value{
  The function \code{PEMGUI} returns an instance of the
  \code{\linkS4class{bayesmodelClass}} class containing
  following information \item{\code{convergence}}{logical,
  whether the model has converged (assessed by the user)}
  \item{\code{results}}{data frame containing statitsics of
  the posterior distribution} \item{\code{jointpost}}{data
  frame giving the joint posterior probability
  distribution} \item{\code{nodes}}{names of the parameters
  jointly estimated by the Bayes model}
  \item{\code{model}}{model in BRugs/Winbugs syntax as a
  character string} \item{\code{chains}}{number of
  independent MCMC chains} \item{\code{burn}}{length of
  burn-in period} \item{\code{update}}{length of update
  iterations for estimation}
}
\description{
  This function provides a GUI for the function
  \link{rrisk.BayesPEM}.
}
\details{
  The diagnostic parameters \code{se} and \code{sp} are
  defined at the pool level, equivalent to
  \code{missclass='pool'} in \code{rrisk.BayesPEM}
  function. See \code{\link{rrisk.BayesPEM}} for more
  details.
}
\note{
  The convergence of the model is assessed by the user
  using diagnostic plots provided by the \pkg{BRugs}
  package.
}
\examples{
#------------------------------------------
# Example of PEM model. Without parameters,
# the input fields will show default values
#------------------------------------------
\donttest{mod <- PEMGUI()}
}
\seealso{
  \code{\link{rrisk.BayesPEM}}
}
\keyword{manip}

