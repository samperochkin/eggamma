#' Control Parameters for dggamma and nllggamma
#'
#' Creates a list of control parameters for \code{dggamma} and \code{nllggamma}.
#'
#' @param ... Named parameters to override defaults.
#' @return A list with control parameters.
#' @examples
#' expansionControl(epsilon = c(1e-2, 1/3), kappa = 10L)
#' @export
expansionControl <- function(...){
  con <- list(epsilon = c(1e-1, 1/2), kappa = 12L)
  modifyList(con, list(...))
}

#' Generalized Gamma Density Function
#'
#' Computes the density of the generalized gamma distribution using a series expansion.
#'
#' @param x Numeric vector of values where the density is evaluated.
#' @param mu Positive numeric, location parameter.
#' @param sigma Positive numeric, scale parameter.
#' @param nu Numeric, shape parameter.
#' @param return_log Logical, whether to return the log-density. Default is FALSE.
#' @param negate_loglik Logical, whether to negate the log-likelihood. Default is FALSE.
#' @param control A list of control parameters, typically created using \code{expansionControl()}.
#'
#' @details
#' The control parameters are \code{epsilon} and \code{kappa}. \code{epsilon} is a length two numeric vector 
#' indicating the intervals \code{c(-epsilon[2],epsilon[1])} and \code{c(epsilon[1],epsilon[2])} where the
#' transition is made from the standard generalized gamma log-lik to its approximation. \code{kappa} 
#' determines the maximal power of \code{nu} to consider in the Taylor and Stirling expansions formulas.

#' 
#' @return A numeric vector of density values.
#' @examples
#' dggamma(x = c(1, 2, 3), mu = 2, sigma = 1, nu = 0.5)
#' dggamma(x = 1:5, mu = 3, sigma = 2, nu = -0.2, return_log = TRUE)
#' @export
dggamma <- function(x, mu, sigma, nu, 
                return_log = FALSE, negate_loglik = FALSE, 
                control = expansionControl()) {
  
  # Basic conditions on parameters
  if(any(x <= 0)) stop("x must be positive.")
  if(any(mu <= 0)) stop("mu must be positive.")
  if(any(sigma <= 0)) stop("sigma must be positive.")
  
  # Conditions on parameters' length
  lens <- sapply(list(x,mu,sigma,nu), length)
  n <- max(lens)
  if(!all(lens %in% c(1,n))) stop("the lengths of x, mu, sigma and nu don't match (length one is allowed).")
  
  # We always want that a vector of length n
  if(lens[1] == 1) x <- rep(x, n)
  if(lens[2] == 1) mu <- rep(mu, n)
  
  # We want both these vectors of length n if any one of the two is
  if(any(lens[3:4] > 1)){
    if(lens[3] == 1) sigma <- rep(sigma, n)
    if(lens[4] == 1) nu <- rep(nu, n)
  }
  
  list2env(control, envir = environment())
  if(!exists("kappa")) stop("No 'kappa' provided (in control), but required for method='expansion'.")
  if(!exists("epsilon")) stop("No 'epsilon' provided (in control), but required for method='expansion'.")
  if(length(epsilon) < 2) stop("'epsilon' provided (in control) must be of length 2 for method='expansion'.")
  epsilon <- epsilon[1:2]
  
  return(dggamma_exp_cpp(x, mu, sigma, nu, epsilon, kappa, return_log, negate_loglik))
}

#' Negative Log-Likelihood for Generalized Gamma
#'
#' Computes the negative log-likelihood (w/ gradient and Hessian) for the generalized gamma distribution.
#'
#' @param p A numeric vector of parameters.
#' @param x Numeric vector of data.
#' @param links Character vector (either "identity", "log", or "logit") specifying the three link functions for parameters.
#' @param control A list of control parameters, typically created using \code{expansionControl()}.
#'
#' @return The negative log-likelihood, with attributes for the gradient and Hessian matrix.
#' 
#' @details
#' This function is meant to be used with nlm to fit an extended generalized distribution.
#' 
#' The parameters \code{p} are those obtained by applying the link functions to mu, sigma, and nu.
#' In other words, using \code{links = rep("identity", 3)} implies that \code{p = c(mu, sigma, nu)}.
#' Default is \code{links = c("log", "log", "identity")}.
#' 
#' The control parameters are \code{epsilon} and \code{kappa}. \code{epsilon} is a length two numeric vector 
#' indicating the intervals \code{c(-epsilon[2],epsilon[1])} and \code{c(epsilon[1],epsilon[2])} where the
#' transition is made from the standard generalized gamma negative log-lik to its approximation. \code{kappa} 
#' determines the maximal power of \code{nu} to consider in the Taylor and Stirling expansions formulas.
#' 
#' @examples
#' nllggamma(p = c(1, 0, 0), x = c(1, 2, 3), links = c("identity", "log", "logit"))
#' @useDynLib eggamma
#' @export
nllggamma <- function(p, x, links, w = 1, control = expansionControl()) {
  
  # links_fun <- list(\(x) x, log, \(x) 1/(1+exp(-x)))
  ilinks_fun <- list(\(x) x, exp, \(x) log(x/(1-x)))
  links_id <- match(links, c("identity", "log", "logit"))
  p <- sapply(1:3, \(i) ilinks_fun[[links_id[i]]](p[i]))
  mu <- p[1]; sigma <- p[2]; nu <- p[3]
  
  # Basic conditions on parameters
  if(any(x <= 0)) stop("x must be positive.")
  if(any(mu <= 0)) stop("mu must be positive.")
  if(any(sigma <= 0)) stop("sigma must be positive.")
  
  # Conditions on parameters' length
  lens <- sapply(list(x,mu,sigma,nu,w), length)
  n <- max(lens)
  if(!all(lens %in% c(1,n))) stop("the lengths of x, mu, sigma and nu don't match (length one is allowed).")
  
  # We always want that a vector of length n
  if(lens[1] == 1) x <- rep(x, n)
  if(lens[2] == 1) mu <- rep(mu, n)
  if(lens[3] == 1) sigma <- rep(sigma, n)
  if(lens[4] == 1) nu <- rep(nu, n)
  if(lens[5] == 1) w <- rep(w, n)
  
  list2env(control, envir = environment())
  
  h1g_fun <- list(\(x) 1, \(x) x, \(x) exp(-x)/(1+exp(-x))^2)
  h2g_fun <- list(\(x) 0, \(x) x, \(x) exp(-2*x)*(exp(x)-1)*(exp(x)+1)^3)
  h1g <- sapply(1:3, \(i) h1g_fun[[links_id[i]]](p[i]))
  h2g <- sapply(1:3, \(i) h2g_fun[[links_id[i]]](p[i]))
  
  if(!exists("kappa")) stop("No 'kappa' provided (in control), but required for method='expansion'.")
  if(!exists("epsilon")) stop("No 'epsilon' provided (in control), but required for method='expansion'.")
  if(length(epsilon) < 2) stop("'epsilon' provided (in control) must be of length 2 for method='expansion'.")
  epsilon <- epsilon[1:2]
  
  nll_plus <- nllggamma_cpp(x, mu, sigma, nu, w, epsilon, kappa)
  nll <- nll_plus[1]
  attr(nll, "gradient") <- nll_plus[2:4] * h1g
  H <- diag(nll_plus[5:7])
  H[cbind(c(1,1,2,2,3,3),c(2,3,3,1,1,2))] <- rep(nll_plus[8:10], times=2)
  H <- H * tcrossprod(h1g) + diag(nll_plus[2:4]*h2g)
  attr(nll, "hessian") <- H
  return(nll)
}



#' Generalized Gamma Quantile Function
#'
#' Computes the quantiles of the generalized gamma distribution using a simple linear interpolation around nu=0.
#'
#' @param p vector of probabilities.
#' @param mu Positive numeric, location parameter.
#' @param sigma Positive numeric, scale parameter.
#' @param nu Numeric, shape parameter.
#' @param log.p logical; if \code{TRUE}, probabilities \code{p} are returned as \code{log(p)}.
#' @param lower.tail Logical; defined as usual.
#' @param epsilon Positive number defining when to start using linear interpolation (from \code{abs(nu)=epsilon} to \code{nu=0}).
#' 
#' @return A numeric vector of quantiles.
#' @export
qggamma <- function (p, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, 
          log.p = FALSE, epsilon = 1e-4) 
{
  if (any(mu < 0)) stop("mu must be positive.")
  if (any(sigma < 0)) stop("sigma must be positive.")
  if (log.p == TRUE) p <- exp(p)
  if (!lower.tail) p <- 1 - p
  if (any(p < 0) | any(p > 1)) stop("p must be between 0 and 1.")
  
  # could certainly be optimized 
  # Conditions on parameters' length
  lens <- sapply(list(p,mu,sigma,nu), length)
  n <- max(lens)
  if(!all(lens %in% c(1,n))) stop("The lengths of x, mu, sigma and nu don't match (length one is allowed).")
  if(lens[1] == 1) p <- rep(p, n)
  if(lens[2] == 1) mu <- rep(mu, n)
  if(lens[3] == 1) sigma <- rep(sigma, n)
  if(lens[4] == 1) nu <- rep(nu, n)
  
  p <- ifelse(nu > 0, p, 1 - p)
  w <- pmax(pmin(abs(nu)/epsilon,1),0)
  nu <- sign(nu)*pmax(epsilon, abs(nu))
  nu[nu == 0] <- epsilon
  xi_inv <- (sigma*nu)^2
  q1 <- mu*qgamma(p, shape=1/xi_inv, scale=xi_inv)^(1/nu) 
  q2 <- qlnorm(p, log(mu), sigma)
  
  w*q1 + (1-w)*q2
}




#' #' Generalized Gamma CDF
#' #'
#' #' Computes the cdf of the generalized gamma distribution using a simple linear interpolation around nu=0.
#' #'
#' #' @param x vector of values.
#' #' @param mu Positive numeric, location parameter.
#' #' @param sigma Positive numeric, scale parameter.
#' #' @param nu Numeric, shape parameter.
#' #' @param log.p logical; if \code{TRUE}, probabilities \code{p} are returned as \code{log(p)}.
#' #' @param lower.tail Logical; defined as usual.
#' #' @param epsilon Positive number defining when to start using linear interpolation (from \code{abs(nu)=epsilon} to \code{nu=0}).
#' #' 
#' #' @return A numeric vector of quantiles.
#' #' @export
#' pggamma <- function (x, mu = 1, sigma = 0.5, nu = 1, lower.tail = TRUE, 
#'                      log.p = FALSE, epsilon = 1e-4) 
#' {
#'   if (any(mu < 0)) stop("mu must be positive.")
#'   if (any(sigma < 0)) stop("sigma must be positive.")
#'   if (log.p == TRUE) p <- exp(p)
#'   if (!lower.tail) p <- 1 - p
#'   if (any(p < 0) | any(p > 1)) stop("p must be between 0 and 1.")
#'   
#'   # could certainly be optimized 
#'   # Conditions on parameters' length
#'   lens <- sapply(list(x,mu,sigma,nu), length)
#'   n <- max(lens)
#'   if(!all(lens %in% c(1,n))) stop("The lengths of x, mu, sigma and nu don't match (length one is allowed).")
#'   if(lens[1] == 1) x <- rep(x, n)
#'   if(lens[2] == 1) mu <- rep(mu, n)
#'   if(lens[3] == 1) sigma <- rep(sigma, n)
#'   if(lens[4] == 1) nu <- rep(nu, n)
#'   
#'   w <- pmax(pmin(abs(nu)/epsilon,1),0)
#'   nu <- sign(nu)*pmax(epsilon, abs(nu))
#'   nu[nu == 0] <- epsilon
#'   xi_inv <- (sigma*nu)^2
#'   ## ************ q1 <- pGA(z, mu = 1, sigma = sigma * abs(nu), lower.tail = (nu < 0) - lower.tail, log.p = log.p)
#'   q2 <- qlnorm(p, log(mu), sigma)
#'   
#'   w*q1 + (1-w)*q2
#' }
