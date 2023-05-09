\name{eb}
\alias{eb}
\title{
Function for Entropy Balancing
}
\description{
This function is called internally by \code{\link{ebalance}} and \code{\link{ebalance.trim}} to implement entropy balancing. This function would normally not be called manually by a user.
}
\usage{
eb(tr.total = tr.total, co.x = co.x,
   coefs = coefs, base.weight = base.weight, 
   max.iterations = max.iterations, 
   constraint.tolerance = constraint.tolerance, 
   print.level = print.level)
}

\arguments{
  \item{tr.total}{
NA
}
  \item{co.x}{
NA
}
  \item{coefs}{
NA
}
  \item{base.weight}{
NA
}
  \item{max.iterations}{
NA
}
  \item{constraint.tolerance}{
NA
}
  \item{print.level}{
NA
}
}

\value{ A list containing the results from the algorithm.
}

\author{
Jens Hainmueller
}

\seealso{
\code{ebalance}, \code{ebalance.trim}
}
\examples{
##---- NA -----
}
