\name{line.searcher}
\alias{line.searcher}
\title{
Optimal step length search for entropy balancing algorithm 
}
\description{
Function called internally by \code{ebalance} and \code{ebalance.trim} to compute optimal step length for entropy balancing algorithm. This function would normally not be called manually by a user.  
}
\usage{
line.searcher(Base.weight, Co.x, 
Tr.total, coefs, Newton, ss)
}
\arguments{
  \item{Base.weight}{
NA
}
  \item{Co.x}{
NA
}
  \item{Tr.total}{
NA
}
  \item{coefs}{
NA
}
  \item{Newton}{
NA
}
  \item{ss}{
NA
}
}

\value{
A list with the results from the search.
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