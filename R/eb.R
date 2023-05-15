# %% ####################################################
#' internal function to compute entropy balancing weights
#' @param tr.total treatment moments
#' @param co.x covariate matrix
#' @param coefs starting coefs
#' @param base.weight base weight
#' @param max.iterations max number of iterations
#' @param constraint.tolerance tolerance for constraint
#' @param print.level print level for debugging
#' @return list with results
eb = function(tr.total = tr.total,
               co.x = co.x,
               coefs = coefs,
               base.weight = base.weight,
               max.iterations = max.iterations,
               constraint.tolerance = constraint.tolerance,
               print.level = print.level) {
  converged = FALSE
  for (iter in 1:max.iterations) {
    weights.temp = c(exp(co.x %*% coefs))
    weights.ebal = weights.temp * base.weight
    co.x.agg = c(weights.ebal %*% co.x)
    gradient = co.x.agg - tr.total
    if (max(abs(gradient)) < constraint.tolerance) {
      converged = TRUE
      break
    }
    if (print.level >= 2) {
      cat("Iteration", iter, "maximum deviation is =", format(max(abs(gradient)), digits = 4), "\n")
    }
    hessian = t(co.x) %*% (weights.ebal * co.x)
    Coefs = coefs
    newton = solve(hessian, gradient)
    coefs = coefs - newton
    loss.new = line.searcher(Base.weight = base.weight, Co.x = co.x, Tr.total = tr.total, coefs = Coefs, Newton = newton, ss = 1)
    loss.old = line.searcher(Base.weight = base.weight, Co.x = co.x, Tr.total = tr.total, coefs = Coefs, Newton = newton, ss = 0)
    if (print.level >= 3) {
      cat("new loss", loss.new, "old loss=", loss.old, "\n")
    }
    if (loss.old <= loss.new) {
      ss.out = optimize(line.searcher,
        lower = .00001, upper = 1, maximum = FALSE,
        Base.weight = base.weight, Co.x = co.x, Tr.total = tr.total, coefs = Coefs, Newton = newton
      )

      if (print.level >= 3) {
        cat("LS Step Length is ", ss.out$minimum, "\n")
      }
      if (print.level >= 3) {
        cat("Loss is", ss.out$objective, "\n")
      }
      coefs = Coefs - ss.out$minimum * solve(hessian, gradient)
    }
  }
  if (print.level >= 1 && converged) {
    cat("Converged within tolerance \n")
  }
  return(
    list(
      maxdiff = max(abs(gradient)),
      coefs = coefs,
      Weights.ebal = weights.ebal,
      converged = converged
    )
  )
}

# function to conduct line search for optimal step length
line.searcher = function(Base.weight,
                          Co.x,
                          Tr.total,
                          coefs,
                          Newton,
                          ss) {
  weights.temp = c(exp(Co.x %*% (coefs - (ss * Newton))))
  # weights.temp[is.infinite(weights.temp)] = 100
  weights.temp = weights.temp * Base.weight
  Co.x.agg = c(weights.temp %*% Co.x)
  maxdiff = max(abs(Co.x.agg - Tr.total))
  return(maxdiff)
}
