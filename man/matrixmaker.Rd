\name{matrixmaker}
\alias{matrixmaker}


\title{
Generate Matrix of One-way Interactions and Squared Terms

}
\description{
Takes a matrix of covariates and generates a new matrix that contains the original covariates, all one-way interaction terms, and all squared terms.
}
\usage{
matrixmaker(mat)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{mat}{n by k numeric matrix of covariates.
}
}

\value{
n by (k*(k+1))/2 +1) matrix of covariates with original covariates, all one-way interaction terms, and all squared terms. 
}

\author{
Jens Hainmueller
}

\seealso{
See \code{\link{getsquares}}
}

\examples{

# create toy matrix
mold <- replicate(3,rnorm(50))
colnames(mold) <- paste("x",1:3,sep="")
head(mold)
# create new matrix
mnew <- matrixmaker(mold)
head(mnew)

}