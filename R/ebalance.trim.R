# %% ####################################################
#' Trimming Weights for Entropy Balancing
#' @description
#' Function to trim weights obtained from entropy balancing. It takes
#' the output from a call to ebalance and trims the data weights
#' (subject to the moment conditions) such that the ratio of the
#' maximum/minimum weight to the mean weight is reduced to satisfy a
#' user specified target. If no user target is specified the maximum
#' weight ratio is automatically trimmed as far as is feasible given the
#' data.
#'
#' @param ebalanceobj	An object from a call to ebalance.
#' @param max.weightOptional target for the ratio of the maximum to mean weight.
#' @param min.weight Optional target for the ratio of the minimum to mean weight.
#' @param max.trim.iterations Maximum number of trimming iterations.
#' @param max.weight.increment Increment for iterative trimming of the ratio of the maximum to mean weight (a scalar between 0-1, .92 indicates that the attempted reduction in the max ratio is 8 percent).
#' @param min.weight.increment Increment for iterative trimming of the ratio of the minimum to mean weight (a scalar > 1, 1.08 indicates that the attempted reduction in the max ratio is 8 percent).
#' @param print.level Controls the level of printing: 0 (normal printing), 2 (detailed), and 3 (very detailed).

#' @return A list object of class ebalance.trim with the following elements:
#' \item{target.margins} A vector that contains the target moments coded from the covariate distributions of the treatment group.
#' \item{co.xdata} A matrix that contains the covariate data from the control group.
#' \item{w} A vector that contains the control group weights assigned by trimming entropy balancing algorithm.
#' \item{coefs} A vector that contains coefficients from the reweighting algorithm.
#' \item{maxdiff} A scalar that contains the maximum deviation between the moments of the reweighted data and the target moments.
#' \item{norm.constant} Normalizing constant used.
#' \item{constraint.tolerance} The tolerance level used for the balance constraints.
#' \item{max.iterations} Maximum number of trimming iterations used.
#' \item{base.weight} The base weight used.
#' \item{converged} Logical flag if algorithm converged within tolerance.
#' @examples
#' # create toy data: treatment indicator and three covariates X1-3
#' treatment   = c(rep(0,50),rep(1,30))
#' X           = rbind(replicate(3,rnorm(50,0)),replicate(3,rnorm(30,.5)))
#' colnames(X) = paste("x",1:3,sep="")
#' # entropy balancing
#' eb.out = ebalance(Treatment=treatment, X=X)
#' # means in treatment group data
#' apply(X[treatment==1,],2,mean)
#' # means in reweighted control group data
#' apply(X[treatment==0,],2,weighted.mean,w=eb.out$w)
#' # means in raw data control group data
#' apply(X[treatment==0,],2,mean)
#' eb.out.tr = ebalance.trim(eb.out)
#' eb.out.tr %>% str
#' # means in reweighted control group data
#' apply(X[treatment==0,],2,weighted.mean,w=eb.out.tr$w)
#' # untrimmed and trimmed weights
#' round(summary(eb.out$w),2)
#' round(summary(eb.out.tr$w),2)
#' @references Hainmueller, J. (2012) 'Entropy Balancing for Causal Effects: A Multivariate Reweighting Method to Produce Balanced Samples in Observational Studies', Political Analysis (Winter 2012) 20 (1): 25–46.
#' @references Zaslavsky, A. (1988), 'Representing local reweighting area adjustments by of households', Survey Methodology 14(2), 265–288.
#' @references Ireland, C. and Kullback, S. (1968), 'Contingency tables with given marginals', Biometrika 55, 179–188.
#' @references Kullback, S. (1959), Information Theory and Statistics, Wiley, NY.
#' @export

ebalance.trim = function(ebalanceobj,
                         max.weight = NULL,
                         min.weight = 0,
                         max.trim.iterations = 200,
                         max.weight.increment = .92,
                         min.weight.increment = 1.08,
                         print.level = 0) {
  if (is(ebalanceobj, "ebalance") == FALSE) {
    stop("ebalanceobj must be an ebalance object from a call to ebalance()")
  }
  minimization = FALSE
  if (is.null(max.weight)) {
    max.weight = max(ebalanceobj$w / mean(ebalanceobj$w))
    minimization = TRUE
  }
  if (length(max.weight) != 1) stop("length(max.weight) != 1")
  if (length(min.weight) != 1) stop("length(min.weight) != 1")
  # starting setup for trimming
  w.trimming = 1
  coefs = ebalanceobj$coefs
  ### Trimming to user supplied max weight or trim once starting from max weight obtained from distribution)
  for (iter.trim in 1:max.trim.iterations) {
    if (minimization == FALSE) {
      cat("Trim iteration", format(iter.trim, digits = 3), "\n")
    }
    eb.out = eb(
      tr.total = ebalanceobj$target.margins,
      co.x = ebalanceobj$co.xdata,
      coefs = coefs,
      base.weight = ebalanceobj$w * w.trimming,
      max.iterations = ebalanceobj$max.iterations,
      constraint.tolerance = ebalanceobj$constraint.tolerance,
      print.level = print.level
    )
    weights.ratio = eb.out$Weights.ebal / mean(eb.out$Weights.ebal)
    coefs = eb.out$coefs
    if (max(weights.ratio) <= max.weight && min(weights.ratio) >= min.weight) {
      if (minimization == FALSE) {
        cat("Converged within tolerance \n")
      }
      break
    }
    w.trimming = w.trimming * ifelse(weights.ratio > max.weight, w.trimming * ((max.weight * max.weight.increment) / weights.ratio), 1)
    if (min.weight) w.trimming = w.trimming *
      ifelse(weights.ratio < min.weight, w.trimming * ((min.weight * min.weight.increment) / weights.ratio), 1)
  }

  ### automated trimming to minimize max weight
  if (minimization == TRUE) {
    cat("Automated trimmig of max weight ratio \n")
    for (iter.max.weight in 1:max.trim.iterations) {
      max.weight.old = max.weight
      weights.ratio.old = max(eb.out$Weights.ebal / mean(eb.out$Weights.ebal))
      eb.out.old = eb.out # store old weights
      if (print.level >= 0) {
        cat("Trim iteration", format(iter.max.weight, digits = 3), "Max Weight Ratio:", format(weights.ratio.old, digits = 4), "\n")
      }
      IsError = try(
        {
          for (iter.trim in 1:max.trim.iterations) {
            eb.out = eb(
              tr.total = ebalanceobj$target.margins,
              co.x = ebalanceobj$co.xdata,
              coefs = coefs,
              base.weight = ebalanceobj$w * w.trimming,
              max.iterations = ebalanceobj$max.iterations,
              constraint.tolerance = ebalanceobj$constraint.tolerance,
              print.level = 0
            )
            weights.ratio = eb.out$Weights.ebal / mean(eb.out$Weights.ebal)
            coefs = eb.out$coefs
            if (max(weights.ratio) <= max.weight && min(weights.ratio) >= min.weight) break
            w.trimming = w.trimming * ifelse(weights.ratio > max.weight, w.trimming * ((max.weight * max.weight.increment) / weights.ratio), 1)
            if (min.weight) w.trimming = w.trimming *
              ifelse(weights.ratio < min.weight, w.trimming * ((min.weight * min.weight.increment) / weights.ratio), 1)
          }
        }, silent = TRUE
      ) # end try call

      if (is(IsError, "try-error")) {
        if (print.level >= 2) {
          cat("no further decrease in max weight ratio \n")
        }
        break
      }
      if (weights.ratio.old < max(weights.ratio)) {
        if (print.level >= 2) {
          cat("no further decrease in max weight ratio \n")
        }
        break
      }
      # exit loop if alog doesn converge
      max.weight = max.weight * max.weight.increment # otherwise lower max weight and try again
      # cat("Moving into next minimisation loop with old max weight ratio:",max(weights.ratio),"and new attempted max weight ratio",max.weight,"\n")
    } # end max weight loop
    # play back the results
    eb.out = eb.out.old
    if (eb.out$converged == TRUE) {
      cat("Converged within tolerance \n")
    }
  } # end automated if

  z = list(
    target.margins = ebalanceobj$tr.total,
    co.xdata = ebalanceobj$co.xdata,
    w = eb.out$Weights.ebal,
    coefs = eb.out$coefs,
    maxdiff = eb.out$maxdiff,
    norm.constant = ebalanceobj$norm.constant,
    constraint.tolerance = ebalanceobj$constraint.tolerance,
    max.iterations = ebalanceobj$max.iterations,
    base.weight = ebalanceobj$base.weight,
    converged = eb.out$converged
  )

  class(z) = "ebalance.trim"
  return(z)
}
