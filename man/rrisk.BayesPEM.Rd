\name{rrisk.BayesPEM}
\alias{rrisk.BayesPEM}
\title{Bayesian Prevalence estimation under misclassification (PEM)}
\usage{
  rrisk.BayesPEM(x, n, k, simulation=FALSE, prior.pi,
    prior.se, prior.sp, misclass="pool",chains=3,
    burn=1000, thin=1, update=10000, workdir=getwd(),
    plots=FALSE)
}
\arguments{
  \item{x}{scalar value for number of pools (\code{k>1}) or
  individual outcomes (\code{k=1}) with positive test
  result}

  \item{n}{scalar value for number of pools tested
  (\code{k>1}) or the sample size in application study
  (\code{k=1})}

  \item{k}{scalar value for number of individual samples
  physically combined into one pool; set \code{k>1} for
  pooled sampling and \code{k=1} for individual sampling}

  \item{simulation}{logical, value \code{TRUE} means the
  function will be called within any simulation routine, in
  this case the graphical diagnostic interface will not be
  invoked (default \code{FALSE})}

  \item{prior.pi}{numeric vector containing parameters of a
  beta distribution as prior for prevalence \code{pi}, e.g.
  \code{pi} \eqn{\sim} \code{prior.pi(*,*)=beta(*,*)}}

  \item{prior.se}{numeric vector containing parameters of a
  beta distribution as prior for sensitivity \code{se},
  e.g. \code{se} \eqn{\sim} \code{prior.se(*,*)=beta(*,*)}}

  \item{prior.sp}{numeric vector containing parameters of a
  beta distribution as prior for specificity \code{sp},
  e.g. \code{sp} \eqn{\sim} \code{prior.sp(*,*)=beta(*,*)}}

  \item{misclass}{character with legal character entries
  \code{pool}, \code{individual} or \code{compare}; ignored
  if k=1}

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
  The function \code{rrisk.BayesPEM} returns an instance of
  the \code{\linkS4class{bayesmodelClass}} class containing
  following informations \item{\code{convergence}}{logical,
  whether the model has converged (assessed by the user)}
  \item{\code{results}}{data frame containing statistics of
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
  Bayesian PEM models provide the posterior distribution
  for the true prevalence (\code{pi}), diagnostic
  sensitivity (\code{se}) and specificity (\code{sp}) for a
  given empirical prevalence estimate using physically
  pooled samples (if \code{k>1}) and priors for the model
  parameters. The misclassification parameters (\code{se}
  and \code{sp}) can be specified at the level of the pool
  or individual level of testing. On the other side, the
  function estimates the true prevalence based on the
  results (\code{x/n}) of an application study with
  individual samples (if \code{k=1}) using a diagnostic
  test, for which some prior information on sensitivity and
  specificity is available.
}
\details{
  The Bayesian model for estimation prevalence, sensitivity
  and specificity has in BRugs/Winbugs syntax following
  form for misclassification at the pool-level (\code{k>1}
  and \code{misclass="pool"}) \preformatted{model{

  pi ~ dbeta(prior.pi[1],prior.pi[2])

  se ~ dbeta(prior.se[1],prior.se[2])

  sp ~ dbeta(prior.sp[1],prior.sp[2])

  p.neg <- pow(1-pi,k)

  p <- (1-p.neg)*se + p.neg*(1-sp)

  x ~ dbin(p,n)

  }} for misclassifications at the individual level
  (\code{k>1} and \code{misclass="individual"})
  \preformatted{model{

  pi ~ dbeta(prior.pi[1],prior.pi[2])

  se ~ dbeta(prior.se[1],prior.se[2])

  sp ~ dbeta(prior.sp[1],prior.sp[2])

  ap <- pi*se + (1-pi)*(1-sp)

  p <- 1- pow(1-ap,k)

  x ~ dbin(p,n)

  }} and for comparison (\code{k>1}) \preformatted{model{

  pi1 ~ dbeta(prior.pi[1],prior.pi[2])

  pi2 ~ dbeta(prior.pi[1],prior.pi[2])

  se ~ dbeta(prior.se[1],prior.se[2])

  sp ~ dbeta(prior.sp[1],prior.sp[2])

  x1 <- x

  x2 <- x

  p.neg <- pow(1-pi1,k)

  p.pos <- (1-p.neg)*se + p.neg*(1-sp)

  x1 ~ dbin(p.pos,n)

  ap <- pi2*se + (1-pi2)*(1-sp)

  p <- 1- pow(1-ap,k)

  x2 ~ dbin(p,n)

  }} The application data (\code{k=1}) has one degree of
  freedom while the underlying model has three unknown
  parameters. Thus, the model is not identifiable and
  informative priors on at least two model parameters are
  required. The Bayesian model for estimation prevalence,
  sensitivity and specificity takes a form
  \preformatted{model{

  x ~ dbin(p,n)

  p <- pi * se + (1-pi) * (1-sp)

  se ~ dbeta(prior.se[1],prior.se[2])

  sp ~ dbeta(prior.sp[1],prior.sp[2])

  pi ~ dbeta(prior.pi[1],prior.pi[2])

  }}
}
\note{
  The convergence of the model is assessed by the user
  using diagnostic plots provided by the \pkg{BRugs}
  package.
}
\examples{
\donttest{
#------------------------------------------
# Example of PEM model (k>1)
#------------------------------------------
pi <- 0.01
se <- 0.96
se.n <- 1000
sp <- 0.99
sp.n <- 1000
n <- sample(10:1000,1,replace=TRUE)  # stochatsic sample size
k <- sample(5:50,1,replace=FALSE)    # stochastic pool size

# Parameters for beta priors
se.a <- se.n*se+1
se.b <- se.n*(1-se)+1
sp.a <- sp.n*sp+1
sp.b <- sp.n*(1-sp)+1

# Random number of positive pools (x) considering uncertainty of se and sp
ap <- pi*se + (1-pi)*(1-sp)
p.pos <- 1-(1-ap)^k
x <- rbinom(1,prob=p.pos,size=n)

# Estimate using Bayes model at individual level
resPEM1 <- rrisk.BayesPEM(x=x, n=n,k=k,
     prior.pi=c(1,1),prior.se=c(se.a,se.b),prior.sp=c(sp.a,sp.b),
     misclass="individual")
resPEM1@results

# Estimate using Bayes model at pool level
resPEM2 <- rrisk.BayesPEM(x=x, n=n,k=k,
     prior.pi=c(1,1),prior.se=c(se.a,se.b),prior.sp=c(sp.a,sp.b),
     misclass="pool")
resPEM2@results

# Estimate using Bayes model compared
resPEM3 <- rrisk.BayesPEM(x=x, n=n,k=k,
     prior.pi=c(1,1),prior.se=c(se.a,se.b),prior.sp=c(sp.a,sp.b),
     misclass="compare")
resPEM3@results

#------------------------------------------
# Example of PEM model (k=1)
#------------------------------------------
# informative priors -> convergence is o.k.
resPEM4<-rrisk.BayesPEM(x=2,n=10,k=1,prior.se=c(12,22),
 prior.sp=c(22,55),prior.pi=c(1,1))
resPEM4@results

# non-informative priors -> convergence is not o.k.
resPEM5<-rrisk.BayesPEM(x=2,n=10,k=1,prior.se=c(1,1),
 prior.sp=c(1,1),prior.pi=c(1,1))
resPEM5@results

# informative priors -> convergence is o.k., without invoking
# graphical diagnostic interface
resPEM6<-rrisk.BayesPEM(x=2,n=10,k=1,prior.se=c(12,22),
 prior.sp=c(22,55),prior.pi=c(1,1))
resPEM6@results
}
}
\references{
  Cowling, D.W., I.A. Gardner and W.O. Johnson (1999).
  Comparison of methods for estimation of individual-level
  prevalence based on pooled samples, Prev.Vet.Med. 39:
  211-225. \cr \cr Rogan, W.J. and B. Gladen (1978).
  Estimating prevalence from the results of a screening
  test. Am. J. Epidemiol. 107: 71-76.
}
\keyword{manip}

