\name{ebalance}
\alias{ebalance}

\title{
Entropy balancing
}
\description{
This function implements entropy balancing, a data preprocessing procedure that allows users to reweight a dataset. The preprocessing is based on a maximum entropy reweighting scheme 
that assigns weights to each unit such that the covariate distributions in the reweighted data satisfy a set of moment conditions specified by the researcher. This can be useful to balance covariate distributions in observational studies with a binary treatment where the control group data can be reweighted to match the covariate moments in the treatment group. Entropy balancing can also be used to reweight a survey sample to known characteristics from a target population. The weights that result from entropy balancing
can be passed to regression or other models to subsequently analyze the reweighted data. 

By default, \code{ebalance} reweights the covariate distributions from a control group to match target moments that are computed from a treatment group such that the reweighted data can be used to analyze the average treatment effect on the treated. 

}
\usage{
ebalance(Treatment, X, base.weight = NULL, 
norm.constant = NULL, coefs = NULL, 
max.iterations = 200, 
constraint.tolerance = 1, print.level = 0)
}

\arguments{
  \item{Treatment}{
A vector indicating the observations which are in the data group that should be reweighted (i.e. the control group) and 
those which are in the data group that should be used to compute the target moments (i.e. the treatment group). 
This can either be a logical vector or a real vector where 0 denotes control observations and 1 denotes treatment observations. 
By default the target moments are computed using the data from all treatment observations. 
}
  \item{X}{
A matrix containing the variables that the researchers wants to include in the reweighting. To adjust the 
means of the covariates, the raw covariates can be included. To adjust the variances of the covariates, squared terms 
of the raw covariates can be included. To adjust co-moments, interaction terms can be included. All columns of 
this matrix must have positive variance and the matrix must be invertible. No missing data is allowed.
}
  \item{base.weight}{
An optional vector of base weights for the maximum entropy reweighting (one weight for each control unit). The default is uniform base weights. 
}
  \item{norm.constant}{
An optional normalizing constant. By default the weights are normalized such that the sum of the weights for the reweighted control group match the number of observations in the treatment grou. 
}
  \item{coefs}{
An optional vector of model coefficients to start the reweighting.
}
  \item{max.iterations}{
Maximum number of iterations that will be run when attempting to reweight the data.
}
  \item{constraint.tolerance}{
This is the tolerance level used by \code{ebalance} to decide if the moments in the reweighted data are equal to the target moments.
}
  \item{print.level}{
Controls the level of printing: 0 (normal printing), 2 (detailed), and 3 (very detailed).
}
}
\details{
If the user wants to pass different target moments these should be entered in a new row in the dataset which contains the target means for the covariate distributions. 
}
\value{
An list object of class \code{ebalance} with the following elements:
  \item{target.margins}{A vector that contains the target moments coded from the covariate distributions of the treatment group.}
  \item{co.xdata}{A matrix that contains the covariate data from the control group.}
  \item{w}{A vector that contains the control group weights assigned by entropy balancing.}
  \item{coefs}{A vector that contains coefficients from the reweighting algorithm.}
  \item{maxdiff}{A scalar that contains the maximum deviation between the moments of the reweighted data and the target moments.}
  \item{constraint.tolerance}{The tolerance level used for the balance constraints.}
  \item{base.weight}{The base weight used.}
  \item{print.level}{The print level used.}
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
}

