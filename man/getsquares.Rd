\name{getsquares}
\alias{getsquares}

\title{
Generate Matrix of Squared Terms
}
\description{
Takes a matrix of covariates and generates a new matrix that contains the original covariates and all squared terms. Squared terms for binary covariates are omitted.
}
\usage{
getsquares(mat)
}
\arguments{
  \item{mat}{n by k numeric matrix of covariates.
  }
}

\value{
n by k*2 numeric matrix that contains the original covariates plus all squared terms.
}

\author{
Jens Hainmueller
}

\seealso{
See \code{\link{matrixmaker}}
}
\examples{
# create toy matrix
mold <- replicate(3,rnorm(50))
colnames(mold) <- paste("x",1:3,sep="")
head(mold)
# create new matrix
mnew <- getsquares(mold)
head(mnew)
}
