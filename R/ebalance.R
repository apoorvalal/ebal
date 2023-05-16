# %% ####################################################
#' Primary Function for Entropy Balancing
#' @description
#' This function implements entropy balancing, a data preprocessing
#' procedure that allows users to reweight a dataset. The preprocessing
#' is based on a maximum entropy reweighting scheme that assigns weights
#' to each unit such that the covariate distributions in the reweighted
#' data satisfy a set of moment conditions specified by the researcher.
#' This can be useful to balance covariate distributions in
#' observational studies with a binary treatment where the control group
#' data can be reweighted to match the covariate moments in the
#' treatment group. Entropy balancing can also be used to reweight a
#' survey sample to known characteristics from a target population. The
#' weights that result from entropy balancing can be passed to
#' regression or other models to subsequently analyze the reweighted
#' data.
#' By default, ebalance reweights the covariate distributions from a
#' control group to match target moments that are computed from a
#' treatment group such that the reweighted data can be used to analyze
#' the average treatment effect on the treated.
#' @param Treatment A vector indicating the observations which are in
#' the data group that should be reweighted (i.e. the control group) and
#' those which are in the data group that should be used to compute the
#' target moments (i.e. the treatment group). This can either be a
#' logical vector or a real vector where 0 denotes control observations
#' and 1 denotes treatment observations. By default the target moments
#' are computed using the data from all treatment observations.
#' @param X A matrix containing the variables that the researchers wants
#' to include in the reweighting. To adjust the means of the covariates,
#' the raw covariates can be included. To adjust the variances of the
#' covariates, squared terms of the raw covariates can be included. To
#' adjust co-moments, interaction terms can be included. All columns of
#' this matrix must have positive variance and the matrix must be
#' invertible. No missing data is allowed.
#' @param base.weight An optional vector of base weights for the
#' maximum entropy reweighting (one weight for each control unit). The
#' default is uniform base weights.
#' @param norm.constant An optional normalizing constant. By default the
#' weights are normalized such that the sum of the weights for the
#' reweighted control group match the number of observations in the
#' treatment group.
#' @param coefs An optional vector of model coefficients to start the
#' reweighting.
#' @param max.iterations Maximum number of iterations that will be run
#' when attempting to reweight the data.
#' @param constraint.tolerance This is the tolerance level used by
#' ebalance to decide if the moments in the reweighted data are equal to
#' the target moments.
#' @param print.level Controls the level of printing: 0 (normal
#' printing), 2 (detailed), and 3 (very detailed).
#' @return An list object of class ebalance with the following elements:
#' \item{target.margins}{A vector that contains the target moments coded from the covariate distributions of the treatment group.}
#' \item{co.xdata}{A matrix that contains the covariate data from the control group.}
#' \item{w}{A vector that contains the control group weights assigned by entropy balancing.}
#' \item{coefs}{A vector that contains coefficients from the reweighting algorithm.}
#' \item{maxdiff}{A scalar that contains the maximum deviation between the moments of the reweighted data and the target moments.}
#' \item{constraint.tolerance}{The tolerance level used for the balance constraints.}
#' \item{base.weight}{The base weight used.}
#' \item{print.level}{The print level used.}
#' \item{converged}{Logical flag if algorithm converged within tolerance.}
#' @examples
#' treatment  =  c(rep(0,50),rep(1,30))
#' X          =  rbind(replicate(3,rnorm(50,0)),replicate(3,rnorm(30,.5)))
#' colnames(X) = paste("x",1:3,sep="")
#' # entropy balancing
#' eb.out  = ebalance(Treatment=treatment, X=X)
#' # means in treatment group data
#' apply(X[treatment==1,],2,mean)
#' # means in reweighted control group data
#' apply(X[treatment==0,],2,weighted.mean,w=eb.out$w)
#' # means in raw data control group data
#' apply(X[treatment==0,],2,mean)
#' @references Hainmueller, J. (2012) 'Entropy Balancing for Causal Effects: A Multivariate Reweighting Method to Produce Balanced Samples in Observational Studies', Political Analysis (Winter 2012) 20 (1): 25–46.
#' @references Zaslavsky, A. (1988), 'Representing local reweighting area adjustments by of households', Survey Methodology 14(2), 265–288.
#' @references Ireland, C. and Kullback, S. (1968), 'Contingency tables with given marginals', Biometrika 55, 179–188.
#' @references Kullback, S. (1959), Information Theory and Statistics, Wiley, NY.
#' @export
ebalance = function(Treatment,
                    X,
                    base.weight = NULL,
                    norm.constant = NULL,
                    coefs = NULL,
                    max.iterations = 200,
                    constraint.tolerance = 1,
                    print.level = 0,
                    method = c("GaussNewton", "AutoDiff")) {
  method = match.arg(method)
  ######################################################################
  # Checks
  ######################################################################
  if (sum(Treatment != 1 & Treatment != 0) > 0) {
    stop("Treatment indicator ('Treatment') must be a logical variable, TRUE (1) or FALSE (0)")
  }
  if (var(Treatment) == 0) {
    stop("Treatment indicator ('Treatment') must contain both treatment and control observations")
  }
  X = as.matrix(X)
  Treatment = as.numeric(Treatment)
  if (sum(is.na(X)) > 0) stop("X contains missing data")
  if (sum(is.na(Treatment)) > 0) stop("Treatment contains missing data")
  if (length(Treatment) != nrow(X)) stop("length(Treatment) != nrow(X)")
  if (length(max.iterations) != 1) stop("length(max.iterations) != 1")
  if (length(constraint.tolerance) != 1) stop("length(constraint.tolerance) != 1")
  # set up elements
  ntreated = sum(Treatment == 1)
  ncontrols = sum(Treatment == 0)
  if (is.null(base.weight)) base.weight = rep(1, ncontrols)
  if (length(base.weight) != ncontrols) {
    stop("length(base.weight) !=  number of controls  sum(Treatment==0)")
  }
  ######################################################################
  # Call optimizer
  ######################################################################
  if (method == "AutoDiff") { # autodiff doesn't need much setup
    target = colMeans(X[Treatment == 1, ])
    source = X[Treatment == 0, ]
    # torch fitter
    eb.out = ebal_torch(
      X0 = source,           # donor moments
      X1 = target,           # target moments
      base_weight = base.weight,
      maxit = max.iterations
    )
    # return list - fewer elements than other option
    z = list(
      target.margins = target,
      co.xdata = source,
      w = eb.out$Weights.ebal,
      coefs = eb.out$coefs,
      max.iterations = max.iterations,
      base.weight = base.weight
    )
    class(z) = "ebalance"
    return(z)
  } else if (method == "GaussNewton") {
    # control units matrix with intercept
    co.x = X[Treatment == 0, ]
    co.x = cbind(rep(1, ncontrols), co.x)
    if (qr(co.x)$rank != ncol(co.x))
      stop("collinearity in covariate matrix for controls (remove collinear covariates)")
    # target moments - sum instead of mean for numerical stability?
    tr.total = apply(as.matrix(X[Treatment == 1, , drop = FALSE]), 2, sum)
    if (is.null(norm.constant)) norm.constant = ntreated
    if (length(norm.constant) != 1) stop("length(norm.constant) != 1")
    tr.total = c(norm.constant, tr.total)
    if (is.null(coefs)) {
      coefs = c(log(tr.total[1] / sum(base.weight)), rep(0, (ncol(co.x) - 1)))
    }
    if (length(coefs) != ncol(co.x)) {
      stop("coefs needs to have same length as number of covariates plus one")
    }
    ## run algo
    eb.out = eb(
      tr.total = tr.total,
      co.x = co.x,
      coefs = coefs,
      base.weight = base.weight,
      max.iterations = max.iterations,
      constraint.tolerance = constraint.tolerance,
      print.level = print.level
    )
    if (eb.out$converged == TRUE & print.level > 0) cat("Converged within tolerance \n")
    z = list(
      target.margins = tr.total,
      co.xdata = co.x,
      w = eb.out$Weights.ebal,
      coefs = eb.out$coefs,
      maxdiff = eb.out$maxdiff,
      norm.constant = norm.constant,
      constraint.tolerance = constraint.tolerance,
      max.iterations = max.iterations,
      base.weight = base.weight,
      print.level = print.level,
      converged = eb.out$converged
    )
    class(z) = "ebalance"
    return(z)
  }
}
