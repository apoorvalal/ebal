\name{ebalance.trim}
\alias{ebalance.trim}
\title{
Trimming of Weights for Entropy Balancing
}
\description{
Function to trim weights obtained from entropy balancing. It takes the output from a call to \code{ebalance} and trims the data weights (subject to the moment conditions) such that the 
ratio of the maximum/minimum weight to the mean weight is 
reduced to satisfy a user specified target. If no user target is specified the maximum weight ratio is automatically trimmed as far as is feasible given the data.
}
\usage{
ebalance.trim(ebalanceobj, max.weight = NULL,
min.weight = 0, max.trim.iterations = 200,
max.weight.increment = 0.92,
min.weight.increment = 1.08,
print.level = 0)
}

\arguments{
  \item{ebalanceobj}{
An object from a call to \code{ebalance}.
}
  \item{max.weight}{
Optional target for the ratio of the maximum to mean weight.
}
  \item{min.weight}{
Optional target for the ratio of the minimum to mean weight.
}         
  \item{max.trim.iterations}{
Maximum number of trimming iterations.
}
  \item{max.weight.increment}{
Increment for iterative trimming of the ratio of the maximum to mean weight (a scalar between 0-1, .92 indicates that 
the attempted reduction in the max ratio is 8 percent). 
}
  \item{min.weight.increment}{
Increment for iterative trimming of the ratio of the minimum to mean weight (a scalar  > 1, 1.08 indicates that 
the attempted reduction in the max ratio is 8 percent).
}
  \item{print.level}{
Controls the level of printing: 0 (normal printing), 2 (detailed), and 3 (very detailed).
}
}

\value{
An list object of class \code{ebalance.trim} with the following elements:
  \item{target.margins}{A vector that contains the target moments coded from the covariate distributions of the treatment group.}
  \item{co.xdata}{A matrix that contains the covariate data from the control group.}
  \item{w}{A vector that contains the control group weights assigned by trimming entropy balancing algorithm.}
  \item{coefs }{A vector that contains coefficients from the reweighting algorithm.}
  \item{maxdiff }{A scalar that contains the maximum deviation between the moments of the reweighted data and the target moments.}
  \item{norm.constant}{Normalizing constant used.}
  \item{constraint.tolerance}{The tolerance level used for the balance constraints.}
  \item{max.iterations}{Maximum number of trimming iterations used.}
  \item{base.weight}{The base weight used.}
  \item{converged}{Logical flag if algorithm converged within tolerance.}
}

\references{

Hainmueller, J. (2012) 'Entropy Balancing for Causal Effects: A Multivariate Reweighting Method to Produce Balanced Samples in Observational Studies', Political Analysis (Winter 2012) 20 (1): 25--46.

Zaslavsky, A. (1988), 'Representing local reweighting area adjustments by of households', Survey Methodology 14(2), 265--288.

Ireland, C. and Kullback, S. (1968), 'Contingency tables with given marginals', Biometrika 55, 179--188.

Kullback, S. (1959), Information Theory and Statistics, Wiley, NY.

}
\author{
Jens Hainmueller
}

\seealso{
Also see \code{\link{ebalance.trim}}.
}

\examples{

# create toy data: treatment indicator and three covariates X1-3
treatment   <- c(rep(0,50),rep(1,30))
X           <- rbind(replicate(3,rnorm(50,0)),replicate(3,rnorm(30,.5)))
colnames(X) <- paste("x",1:3,sep="")

# entropy balancing
eb.out <- ebalance(Treatment=treatment,
                   X=X)
# means in treatment group data
apply(X[treatment==1,],2,mean)
# means in reweighted control group data
apply(X[treatment==0,],2,weighted.mean,w=eb.out$w)
# means in raw data control group data
apply(X[treatment==0,],2,mean)

# trim weights
eb.out.tr <- ebalance.trim(eb.out)
# means in reweighted control group data
apply(X[treatment==0,],2,weighted.mean,w=eb.out.tr$w)

# untrimmed and trimmed weights
round(summary(eb.out$w),2)
round(summary(eb.out.tr$w),2)

}
