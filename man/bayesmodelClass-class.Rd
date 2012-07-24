\name{bayesmodelClass-class}
\Rdversion{1.1}
\docType{class}
\alias{bayesmodelClass-class}

\title{S4 class for displaying output of Bayesian models}
\description{S4 class for displaying the output of the functions \code{\link{rrisk.BayesPEM}}
and \code{\link{rrisk.BayesZIP}}.
}
\section{Objects from the Class}{
Objects can be created by calls of the form \code{new("bayesmodelClass", ...)}.
}
\section{Slots}{
  \describe{
    \item{\code{convergence}:}{Object of class \code{"ANY"}}
    \item{\code{results}:}{Object of class \code{"ANY"}}
    \item{\code{jointpost}:}{Object of class \code{"ANY"}}
    \item{\code{nodes}:}{Object of class \code{"ANY"}}
    \item{\code{model}:}{Object of class \code{"ANY"}}
    \item{\code{chains}:}{Object of class \code{"ANY"}}
    \item{\code{burn}:}{Object of class \code{"ANY"}}
    \item{\code{update}:}{Object of class \code{"ANY"}}
  }
}
\section{Methods}{
  \describe{
    \item{show}{\code{signature(object = "bayesmodelClass")}: ... }
	 }
}
\examples{
new("bayesmodelClass")
new("bayesmodelClass",model="Model info...",nodes="Nodes info...")
showClass("bayesmodelClass")
}
\keyword{classes}
