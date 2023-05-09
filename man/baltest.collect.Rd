\name{baltest.collect}
\alias{baltest.collect}
\title{
Collect Covariate Balance Statistics 
}
\description{
A function that summarizes the covariate balance statistics that are computed by \code{MatchBalance(Matching)} in a balance table.  
}
\usage{
baltest.collect(matchbal.out, var.names, after = TRUE)
}

\arguments{
  \item{matchbal.out}{
An object from a call to \code{MatchBalance(Matching)}
}
  \item{var.names}{
A vector of covariate names.
}
  \item{after}{
A logical flag for whether the results from before or after Matching should be summarized. If TRUE \code{baltest.collect} summarizes the results from the covariate balance checks that \code{MatchBalance} computes in the matched data. If FALSE the results from the balance checks in the unmatched data are used.
}
}
\details{
See \code{MatchBalance(Matching)} for details.
}
\value{
  A matrix that contains the covariate balance statistics in tabular format.
}

\author{
Jens Hainmueller
}


\seealso{
\code{MatchBalance} in the \code{Matching} package.
}
\examples{

## load(Matching) to run this example
## create toy data: one treatment indicator and three covariates X1-3
#dat <- data.frame(treatment=rbinom(50,size=1,prob=.5),replicate(3,rnorm(50)))
#covarsname <- colnames(dat)[-1]

## run balance checks
#mout <- MatchBalance(treatment~X1+X2+X3,data=dat)

## summarize in balance table
#baltest.collect(matchbal.out=mout,var.names=covarsname,after=FALSE)

}
