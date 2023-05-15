# %% ####################################################
#' Generate Matrix of Squared Terms
#' @description Takes a matrix of covariates and generates a new matrix
#' that contains the original covariates and all squared terms. Squared
#' terms for binary covariates are omitted.
#' @param mat	n by k numeric matrix of covariates.
#' @return matrix with additional columns
#' @export
getsquares = function(mat) {
  mat.add = matrix(NA, nrow(mat), 0)
  mat.add.names = c()
  for (i in 1:ncol(mat)) {
    if (dim(table(mat[, i])) > 2) {
      mat.add = cbind(mat.add, mat[, i]^2)
      mat.add.names = c(paste(colnames(mat)[i], ".2", sep = ""), mat.add.names)
    }
  }
  colnames(mat.add) = mat.add.names
  out = cbind(mat, mat.add)
}

# %% ####################################################
#' Generate Matrix of One-way Interactions and Squared Terms
#' @description Takes a matrix of covariates and generates a new matrix
#' that contains the original covariates, all one-way interaction terms,
#' and all squared terms.
#' @param mat	n by k numeric matrix of covariates.
#' @return matrix with additional columns
#' @export
matrixmaker = function(mat) {
  k = ncol(mat)
  obs = nrow(mat)
  out = matrix(NA, obs, ((k * (k + 1)) / 2 + 1))
  count = 0
  for (i in 1:k) {
    for (j in 1:i) {
      count = count + 1
      out[, count] = mat[, i] * mat[, j]
      colnames(out)[c(count, ((k * (k + 1)) / 2 + 1))] =
        c(paste(colnames(mat)[i], ".", colnames(mat)[j], sep = ""), "dummy")
    }
  }
  out = cbind(mat, out[, 1:((k * (k + 1)) / 2)])
  out
}
