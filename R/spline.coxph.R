#' @title Fitting a parametric proportional hazard model with M-spline function
#' @description \code{spline.coxph} fitting a parametric proportional hazard model with M-spline function
#'
#' @importFrom spline.coxph spline.coxph
#' @param t.event observed survival time
#' @param event event Indicator (0 or 1)
#' @param Z covariate vector
#' @param xi1 minimum survival time
#' @param xi3 maximum survival time
#' @param shape select shape of baseline hazard
#' @param p0 initial value of parameters
#' @return results of estimation
#' @importFrom stats nlm
#' @importFrom joint.Cox M.spline
#' @importFrom joint.Cox I.spline
#' @export

spline.coxph <- function (t.event, event, Z, xi1 = min(t.event), xi3 = max(t.event),
                          shape = "bathtub1", p0 = rep(0,1+p))
{
  shape.list = list(
    increase  = c(0.05,0.1,0.15,0.3,0.4),
    constant  = c(0.125,0.25,0.25,0.25,0.125),
    decrease  = c(0.4,0.3,0.15,0.1,0.05),
    unimodal1 = c(0.001,0.014,0.97,0.014,0.001),
    unimodal2 = c(0.001,0.8,0.124,0.074,0.001),
    unimodal3 = c(0.001,0.074,0.124,0.8,0.001),
    bathtub1  = c(0.3,0.1995,0.001,0.1995,0.3),
    bathtub2  = c(0.3,0.001,0.1009,0.299,0.3),
    bathtub3  = c(0.3,0.299,0.1009,0.001,0.3)
  )

  para = shape.list[[shape]]

  d = event
  Z = as.matrix(Z)
  p = ncol(Z)

  l.func = function(phi) {
    b  = phi[2:(1 + p)]
    g  = exp(pmin(phi[1], 500))
    r  = as.vector(M.spline(t.event, xi1 = min(t.event), xi3 = max(t.event))%*%(g * para))
    R  = as.vector(I.spline(t.event, xi1 = min(t.event), xi3 = max(t.event))%*%(g * para))
    bZ = as.numeric(Z %*% b)
    l  = sum(d * (log(r) + bZ))
    l  = l - sum(pmin(exp(bZ) * R, exp(500)))
    -l
  }

  res = nlm(l.func, p = p0, hessian = TRUE)

  beta.est = res$est[2:(1 + p)]
  lam.est  = exp(res$est[1])
  H        = -res$hessian
  V        = solve(-H, tol = 10^(-100))
  beta.se  = sqrt(diag(V)[2:(1 + p)])
  lam.se   = sqrt(lam.est %*% V[1,1] %*% lam.est)
  b.lower  = beta.est - 1.96 * beta.se
  b.upper  = beta.est + 1.96 * beta.se
  l.lower  = lam.est  - 1.96 * lam.se
  l.upper  = lam.est  + 1.96 * lam.se

  beta.res = c(estimate = beta.est, SE = beta.se, Lower = b.lower, Upper = b.upper)
  lam.res  = c(estimate = lam.est , SE = lam.se,  Lower = l.lower, Upper = l.upper)
  list(beta = beta.res, lambda = lam.res)
}

