\name{ZIPGUI}
\alias{ZIPGUI}
\title{Bayes estimation of a zero inflated Poisson (ZIP) model}
\usage{
  ZIPGUI(data=NULL, prior.lambda=c(1,10),
    prior.pi=c(0.8,1), chains=3, burn=1000, update=10000,
    thin=1)
}
\arguments{
  \item{data}{a vector of numeric data, containing zeros,
  and of minimal length 10.}

  \item{prior.lambda}{numeric vector containing minimum and
  maximum of a uniform distribution used as prior for the
  Poisson parameter \code{lambda}, e.g. \cr \code{lambda}
  \eqn{\sim} \code{prior.lambda(*,*)=unif(*,*)}}

  \item{prior.pi}{numeric vector containing parameters of a
  beta distribution describing prior knowledge about
  prevalence (proportion of contaminated samples), e.g. \cr
  \code{pi} \eqn{\sim} \code{prior.pi(*,*)=beta(*,*)}}

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
  The function \code{ZIPGUI} returns an instance of the
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
  \link{rrisk.BayesZIP}.
}
\note{
  The convergence of the model will be entered by the user
  after the simulation process.
}
\examples{
\donttest{
data <- rpois(30, 4)
prior.lambda <- c(1, 10)
prior.pi <- c(0.8, 1)
ZIPGUI(data, prior.lambda, prior.pi)}
}
\references{
  Bohning, D., E. Dietz, P. Schlattman, L. Mendonca, and U.
  Kirchner (1999). The zero-inflated Poisson model and the
  decayed, missing and filled teeth index in dental
  epidemiology. Journal of the Royal Statistical Society,
  Series A 162, 195-209.
}
\seealso{
  \code{\link{rrisk.BayesZIP}}
}
\keyword{manip}

