#' ebal implementation with autodiff via torch
#' @param X0 donor units matrix
#' @param X1 target moments
#' @param base_weight base weight vector
#' @return list with weights and coefficients
#' @import torch
#' @export
ebal_torch = function(X0, X1, base_weight = NULL, maxit = 200) {
  if (is.null(base_weight)) base_weight = rep(1, nrow(X0))
  inp = list(x0 = torch_tensor(X0), x1 = torch_tensor(X1), q = torch_tensor(base_weight))
  # loss fn returns tensor
  ebal_loss_torch = function(lambda) {
    inner_sum = inp$q$dot(torch_exp(-1 * torch_matmul(inp$x0, lambda)))
    loss_value = torch_log(inner_sum) + torch_matmul(lambda, inp$x1)
  }
  # gradient - autograd
  loss_grad = function(lambda) {
    lambda_ad = torch_tensor(lambda, requires_grad = TRUE)
    loss = ebal_loss_torch(lambda_ad)
    grad = autograd_grad(loss, lambda_ad)[[1]]
    as.numeric(grad)
  }
  # call optim - fn needs to have numeric output
  coefs = optim(
    fn = function(x) as.numeric(ebal_loss_torch(x)),
    gr = loss_grad, par = rep(1, ncol(X0)),
    method = "BFGS",
    control = list(maxit = maxit)
  )$par
  # extract weights from lagrangian
  wts = exp(-1 * as.matrix(inp$x0) %*% coefs)
  wts = wts / sum(wts) # normalise
  # return
  return(
    list(
      coefs = coefs,
      Weights.ebal = wts
    )
  )
}
