#' Generic Euler-Maruyama Log-Likelihood for an SDE
#'
#' Computes the approximate log-likelihood of a one-dimensional SDE using
#' the Euler-Maruyama discretization.
#'
#' This function can serve as a baseline for likelihood implementations of a given SDE.
#'
#' @param SDE_object An object of class \code{SDE} containing \code{drift}- and \code{diffusion} functions.
#' @param par Numeric vector. Parameter values to evaluate the likelihood in.
#' @param step_length Numeric vector. Real number that contains the temporal resolution of the trajectory used in the estimation.
#' @param data Numeric vector. Trajectory from the SDE model.

#' @return Numeric. The negative log-likelihood of the observed data under the Euler-Maruyama approximation.
#' @export
EM_likelihood_generic <- function(SDE_object, par, step_length, data) {

  x <- data

  # Validate inputs
  if (!is.numeric(step_length) || length(step_length) != 1 || step_length <= 0) {
    stop("step_length must be a single positive numeric value")
  }

  if (!inherits(SDE_object, "SDE")) {
    stop("SDE_object must be of class 'SDE'")
  }

  if (!is.numeric(x) || length(x) < 2) {
    stop("data must be a numeric vector of length >= 2")
  }

  if (!is.function(SDE_object$drift)) stop("SDE_object$drift must be a function")
  if (!is.function(SDE_object$diffusion)) stop("SDE_object$diffusion must be a function")


  # Euler-Maruyama increments
  X_lower <- x[-length(x)]
  X_upper <- x[-1]

  drift_val <- SDE_object$drift(X_lower, par)
  sigma_val <- SDE_object$diffusion(X_lower, par)

  EM_mean <- X_lower + drift_val * step_length
  EM_sd   <- sigma_val * sqrt(step_length)

  # Return negative log-likelihood
  -sum(stats::dnorm(X_upper, mean = EM_mean, sd = EM_sd, log = TRUE))
}

#' Fit an SDE using Euler-Maruyama Maximum Likelihood
#'
#' Optimizes parameters of an SDE object using the Euler-Maruyama approximation.
#'
#' @param SDE_object An object of class \code{SDE} containing drift and diffusion and potentially data and step_length,
#' @param init_par Numeric vector. Initial guess for the parameters.
#' @param step_length Numeric vector. Optional temporal resolution of the trajectory in the data argument, if the SDE_object for instance does contain a trajectory.
#' @param data Numeric vector. Optional realized trajectory from the model given in the SDE_object, if this does not contain a trajectory, this can be passed explicitly here.
#' @param method Optimization method (default "BFGS").
#' @param ... Additional arguments passed to \code{\link[stats]{optim}}.
#'
#' @return A list with elements `par` (estimated parameters) and `value` (negative log-likelihood).
#' @export
fit_SDE_EM <- function(SDE_object,
                       init_par,
                       step_length = NA,
                       data = NA,
                       method = "BFGS",
                       ...) {

  if (!inherits(SDE_object, "SDE")) stop("SDE_object must be of class 'SDE'")
  if (!is.numeric(SDE_object$X_t) && !is.numeric(data)) stop("Must provide numeric data either via SDE instance or directly with the data-argument.")
  if (!is.numeric(SDE_object$step_length) && !is.numeric(step_length)) stop("Must provide numeric step_length either via SDE instance or directly with the step_length-argument.")
  if (!is.numeric(init_par)) stop("init_par must be numeric")


  # SDE_object$par <- init_par
  if(length(SDE_object$step_length) == 1 && is.na(step_length)){
    step_length <- SDE_object$step_length
  }
  if(length(SDE_object$X_t) > 0 && is.na(data)){
    data <- SDE_object$X_t
  }
  # Wrapper for optim
  # EM_wrapper <- function(par_vec) {
  #   temp_SDE <- SDE_object
  #   temp_SDE$par <- par_vec
  #   EM_likelihood_generic(SDE_object = temp_SDE)
  # }

  res <- stats::optim(
    par = init_par,
    fn = EM_likelihood_generic,
    method = method,
    SDE_object = SDE_object,
    data = data,
    step_length = step_length,
    ...
  )

  return(list(par = res$par,
              value = res$value,
              convergence = res$convergence,
              message = res$message))
}

