\name{rrisk.BayesZIP}
\alias{rrisk.BayesZIP}
\title{Bayes estimation of a zero inflated Poisson (ZIP) model}
\usage{
  rrisk.BayesZIP(data, prior.lambda=c(1,10),
    prior.pi=c(0.8,1), simulation=FALSE, chains=3,
    burn=1000, thin=1, update=10000, workdir=getwd(),
    plots=FALSE)
}
\arguments{
  \item{data}{matrix, data frame or data set with positive
  integers, including zeros and of the minimal length 10}

  \item{prior.lambda}{numeric vector containing minimum and
  maximum of a uniform distribution used as prior for the
  Poisson parameter \code{lambda}, e.g. \cr \code{lambda}
  \eqn{\sim} \code{prior.lambda(*,*)=unif(*,*)}}

  \item{prior.pi}{numeric vector containing parameters of a
  beta distribution describing prior knowledge about
  prevalence (proportion of contaminated samples), e.g. \cr
  \code{pi} \eqn{\sim} \code{prior.pi(*,*)=beta(*,*)}}

  \item{simulation}{logical, value \code{TRUE} means the
  function will be called within any simulation routine, in
  this case the graphical diagnostic interface will not be
  invoked (default \code{FALSE})}

  \item{chains}{positive single numeric value, number of
  independent MCMC chains (default 3)}

  \item{burn}{positive single numeric value, length of the
  burn-in period (default 1000)}

  \item{thin}{positive single numeric value (default 1).
  The samples from every kth iteration will be used for
  inference, where k is the value of thin. Setting
  \code{thin > 1} can help to reduce the autocorrelation in
  the sample.}

  \item{update}{positive single numeric value, length of
  update iterations for estimation (default 10000)}

  \item{workdir}{character string giving working directory
  to store temporary data (default \code{getwd()})}

  \item{plots}{logical, if \code{TRUE} the diagnostic plots
  will be displayed in separate windows}
}
\value{
  The function \code{rrisk.BayesZIP} returns an instance of
  the \code{\linkS4class{bayesmodelClass}} class containing
  following informations \item{\code{convergence}}{logical,
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
  Zero-inflated Poisson data are count data with an excess
  number of zeros. The ZIP model involves the Poisson
  parameter \code{lambda} and the prevalence parameter
  \code{pi}.
}
\details{
  The ZIP model applies to count data and can be
  interpreted as a mixture distribution with one component
  comprising the 'true' zeros and another component of
  Poisson distributed values with density parameter
  \code{lambda}. The prevalence parameter \code{pi} refers
  to the proportion of the second, true non-zero component.
  \cr \cr The Bayesian model for estimation prevalence and
  lambda parameter has in BRugs/Winbugs syntax following
  form \preformatted{model{

  lambda ~ dunif(prior.lambda[1],prior.lambda[2])

  pi ~ dbeta(prior.pi[1],prior.pi[2])

  for (i in 1:n) {

  y[i] ~ dpois(mu[i])

  mu[i] <- I[i] * lambda

  I[i] ~ dbern(pi)

  }

  }}
}
\note{
  The convergence of the model should be checked using the
  diagnostic plots see the package \pkg{BRugs}, see also
  \pkg{zicounts}.
}
\examples{
\donttest{
#------------------------------------------
# Example of ZIP model
#------------------------------------------
# generate ZIP data
pi<-0.01
n<-200
lambda<-3.5
zip.data<-rep(0,n)
zip.data[sample(1:n,n*pi,replace=FALSE)]<-rpois(n*pi,lambda=lambda)

# estimate using Bayes model for zero inflated data
resZIP<-rrisk.BayesZIP(data=zip.data, prior.lambda=c(0,100),prior.pi=c(1,1),
 burn=100,update=1000)
resZIP@results

# estimate using Bayes model for zero inflated data without invoking
# graphical diagnostic interface
rrisk.BayesZIP(data=zip.data, prior.lambda=c(0,100),prior.pi=c(1,1),
 burn=100,update=1000,simulation=TRUE)

# compare with naive results ignoring ZIP model
pi.crude <- sum(zip.data>0)/n
lambda.crude <- mean(zip.data)
print(pi.crude)
print(lambda.crude)
resZIP@results
}
}
\references{
  Bohning, D., E. Dietz, P. Schlattman, L. Mendonca, and U.
  Kirchner (1999). The zero-inflated Poisson model and the
  decayed, missing and filled teeth index in dental
  epidemiology. Journal of the Royal Statistical Society,
  Series A 162, 195-209.
}
\keyword{manip}

