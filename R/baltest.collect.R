# %% ####################################################
#' Collect Covariate Balance Statistics
#' @description A function that summarizes the covariate balance statistics that are computed by MatchBalance(Matching) in a balance table.
#' @param matchbal An object from a call to MatchBalance(Matching)
#' @param var.names A vector of covariate names
#' @param after A logical flag for whether the results from before or after Matching should be summarized. If TRUE baltest.collect summarizes the results from the covariate balance checks that MatchBalance computes in the matched data. If FALSE the results from the balance checks in the unmatched data are used.
#' @return vector of treatment effects
#' @export
baltest.collect = function(matchbal.out, var.names, after = TRUE) {
  storemat = matrix(NA, length(var.names), 10)
  colnames(storemat) = c("mean.Tr", "mean.Co", "sdiff", "sdiff.pooled", "var.ratio", "T pval", "KS pval", "qqmeandiff", "qqmeddiff", "qqmaxdiff")
  rownames(storemat) = var.names

  if (after == FALSE) {
    stopifnot(length(matchbal.out$BeforeMatching) == length(var.names))
    for (i in 1:length(matchbal.out$BeforeMatching)) {
      storemat[i, 1] = matchbal.out$BeforeMatching[[i]]$mean.Tr
      storemat[i, 2] = matchbal.out$BeforeMatching[[i]]$mean.Co
      storemat[i, 3] = matchbal.out$BeforeMatching[[i]]$sdiff
      storemat[i, 4] = matchbal.out$BeforeMatching[[i]]$sdiff.pooled
      storemat[i, 5] = matchbal.out$BeforeMatching[[i]]$var.ratio
      storemat[i, 6] = matchbal.out$BeforeMatching[[i]]$p.value
      TT = try(matchbal.out$BeforeMatching[[i]]$ks$ks.boot.pvalue, silent = TRUE)
      storemat[i, 7] = ifelse(is.null(TT), NA, TT)
      storemat[i, 8] = matchbal.out$BeforeMatching[[i]]$qqsummary$meandiff
      storemat[i, 9] = matchbal.out$BeforeMatching[[i]]$qqsummary$mediandiff
      storemat[i, 10] = matchbal.out$BeforeMatching[[i]]$qqsummary$maxdiff
    }
  } else {
    stopifnot(length(matchbal.out$AfterMatching) == length(var.names))
    for (i in 1:length(matchbal.out$AfterMatching)) {
      storemat[i, 1] = matchbal.out$AfterMatching[[i]]$mean.Tr
      storemat[i, 2] = matchbal.out$AfterMatching[[i]]$mean.Co
      storemat[i, 3] = matchbal.out$AfterMatching[[i]]$sdiff
      storemat[i, 4] = matchbal.out$AfterMatching[[i]]$sdiff.pooled
      storemat[i, 5] = matchbal.out$AfterMatching[[i]]$var.ratio
      storemat[i, 6] = matchbal.out$AfterMatching[[i]]$p.value
      TT = try(matchbal.out$AfterMatching[[i]]$ks$ks.boot.pvalue, silent = TRUE)
      storemat[i, 7] = ifelse(is.null(TT), NA, TT)
      storemat[i, 8] = matchbal.out$AfterMatching[[i]]$qqsummary$meandiff
      storemat[i, 9] = matchbal.out$AfterMatching[[i]]$qqsummary$mediandiff
      storemat[i, 10] = matchbal.out$AfterMatching[[i]]$qqsummary$maxdiff
    }
  }
  return(storemat)
}
