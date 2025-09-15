#' Simulate a Stochastic Saddle-Node Bifurcation
#'
#' This function simulates the trajectory of a stochastic system undergoing a saddle-node bifurcation
#' with time-dependent bifurcation parameter. Different types of stochastic noise can be used, including
#' additive, multiplicative, and heavy-tailed distributions.
#'
#' @param step_length Numeric. The time increment for each simulation step.
#' @param par Numeric vector. Model parameters:
#'   \describe{
#'     \item{A}{Scaling of the process.}
#'     \item{m}{Shift of the process.}
#'     \item{lambda_0}{Initial value of the bifurcation parameter.}
#'     \item{sigma}{Diffusion parameter.}
#'     \item{nu}{(Optional) Exponent controlling the time-dependence of lambda. Defaults to 1.}
#'   }
#' @param X_0 Numeric. Initial value of the state variable. If `NA`, the function sets it to the stable fixed point. Default is `NA`.
#' @param tau Numeric. Time scale over which the bifurcation parameter changes. Default is 100.
#' @param t_0 Numeric. Initial time offset before the bifurcation starts. Default is 10.
#' @param beyond_tipping Numeric. Additional simulation time beyond the bifurcation point. Default is 0.
#' @param sample_method Character. Method used for generating random noise. Passed to `fast_rnorm()`. Default is `"auto"` other valid arguments are: `"stats"` and `"dqrng"`.
#' @param model Character. Type of stochastic noise to use. Options are:
#'   \itemize{
#'     \item `"additive"`: additive Gaussian noise
#'     \item `"sqrt"`: multiplicative noise proportional to sqrt(X)
#'     \item `"linear"`: multiplicative noise proportional to X
#'     \item `"t-dist"`: multiplicative noise proportional to sqrt(X^2 + 1)
#'     \item `"F-dist"`: multiplicative noise proportional to sqrt(X * (X + 1))
#'     \item `"Jacobi"`: multiplicative noise proportional to sqrt(X * (1 - X))
#'   }
#'
#' @return A `data.frame` with the following columns:
#'   \describe{
#'     \item{t}{Time points of the simulation.}
#'     \item{X_t}{Simulated trajectory of the state variable.}
#'     \item{lambda_t}{Time-dependent bifurcation parameter.}
#'     \item{alpha_t}{Central parameter in the estimation of the process.}
#'     \item{mu_t}{Location of the stable fixed point.}
#'   }
#'
#' @details
#' The function discretizes the stochastic differential equation describing a saddle-node bifurcation
#' and applies the chosen stochastic noise term at each time step. The time-dependent bifurcation
#' parameter \code{lambda_t} decreases after \code{t_0} according to a power law with exponent \code{nu}.
#'
#' @examples
#' # Simulate with additive noise
#' simulate_stochastic_saddlenode_bifurcation(
#'   step_length = 0.1,
#'   par = c(A = -1, m = 0, lambda_0 = 1, sigma = 0.1),
#'   tau = 50,
#'   t_0 = 5
#' )
#'
#' @export

simulate_stochastic_saddlenode_bifurcation <- function(step_length,
                                                  par,
                                                  X_0 = NA,
                                                  tau = 100,
                                                  t_0 = 10,
                                                  beyond_tipping = 0,
                                                  sample_method = "auto",
                                                  model = "additive"){
  A        <- par[1]
  m        <- par[2]
  lambda_0 <- par[3]
  sigma    <- par[4]
  nu       <- if(length(par) == 5) par[5] else 1

  # Initial point in fixed point of process if none specified
  if(is.na(X_0)){
    X_0 <- m + ifelse(A >=0, 1, -1) * sqrt(abs(lambda_0 / A))
  }

  #Tipping point with added time if specified
  total_time <- tau + t_0 + beyond_tipping
  N          <- as.integer((total_time) / step_length)

  # Initialize the process

  dW            <- fast_rnorm(N, mean = 0, standard_deviation = sqrt(step_length), method = sample_method)
  time          <- step_length * 0:N
  lambda_t      <- lambda_0 * (1 - (time >= t_0) * (time - t_0) / tau)^nu
  alpha_t       <- 2 * sqrt(abs(A * lambda_t))
  mu_t          <- m + sqrt(abs(lambda_t / A))

  switch(model,
         "additive" = X_t <- updatestep_saddlenode_additive_model(
                                X_0 = X_0,
                                lambda_t = lambda_t,
                                A = A,
                                m = m,
                                sigma = sigma,
                                step_length = step_length,
                                dW = dW),
         "sqrt" = X_t <- updatestep_saddlenode_sqrt_model(
                                X_0 = X_0,
                                lambda_t = lambda_t,
                                A = A,
                                m = m,
                                sigma = sigma,
                                step_length = step_length,
                                dW = dW),
         "linear" = X_t <- updatestep_saddlenode_linear_model(
                                X_0 = X_0,
                                lambda_t = lambda_t,
                                A = A,
                                m = m,
                                sigma = sigma,
                                step_length = step_length,
                                dW = dW),
         "t-dist" = X_t <- updatestep_saddlenode_t_dist(
                           X_0 = X_0,
                           lambda_t = lambda_t,
                           A = A,
                           m = m,
                           sigma = sigma,
                           step_length = step_length,
                           dW = dW),
         "F-dist" = X_t <- updatestep_saddlenode_F_dist(
                           X_0 = X_0,
                           lambda_t = lambda_t,
                           A = A,
                           m = m,
                           sigma = sigma,
                           step_length = step_length,
                           dW = dW),
         "Jacobi" = X_t <- updatestep_saddlenode_jacobi(
                           X_0 = X_0,
                           lambda_t = lambda_t,
                           A = A,
                           m = m,
                           sigma = sigma,
                           step_length = step_length,
                           dW = dW)
         )

  data.frame(t = time,
             X_t,
             lambda_t,
             alpha_t,
             mu_t)
}

#' Simulate a Pearson Diffusion Process
#'
#' This function simulates the trajectory of an Ergodic Pearson diffusion process.
#' This class of diffusion constitutes a highly tractable class of diffusions.
#' Chose between the type of diffusion via the `model` argument.
#' The options are The Ornstein-Uhlenbeck Process, The Square-root Process,
#' The Mean-Reverting Geometric Brownian Motion, The t-diffusion, The F-diffusion,
#' and The Jacobi diffusion.
#'
#' @param step_length Numeric. The time increment for each simulation step.
#' @param par Numeric vector of model parameters:
#'   \describe{
#'     \item{beta}{Drift rate - Strengh of mean-reversion}
#'     \item{mu}{Long-term mean of the process.}
#'     \item{sigma}{Diffusion coefficient.}
#'   }
#' @param X_0 Numeric. Initial value of the process. If `NA`, the function sets it to `mu`. Default is `NA`.
#' @param total_time Numeric. Total simulation time.
#' @param sample_method Character. Method used for generating the Wiener increments. Passed to `fast_rnorm()`.
#'   Default is `"auto"`. Other valid options are `"stats"` and `"dqrng"`.
#' @param model Character. The type of Pearson diffusion to simulate. Valid options are:
#'   \itemize{
#'     \item `"OU"`: Ornstein-Uhlenbeck process (additive Gaussian noise)
#'     \item `"sqrt"`: Square-root process (proportional to sqrt(X)) or CIR-MODEL
#'     \item `"linear"`: Linear noise (noise proportional to X) or GARCH MODEL
#'     \item `"t-dist"`: Multiplicative (noise proportional to sqrt(X^2 + 1)))
#'     \item `"F-dist"`: Multiplicative (noise proportional to sqrt(X * (X + 1)))
#'     \item `"Jacobi"`: Jacobi-diffusion - Multiplicative noise proportional to sqrt(X * (1 - X))
#'   }
#'
#' @return A `data.frame` with the following columns:
#'   \describe{
#'     \item{t}{Time points of the simulation.}
#'     \item{X_t}{Simulated trajectory of the process.}
#'   }
#'
#' @details
#' The function discretizes the chosen Pearson diffusion process and applies the corresponding
#' update rule at each time step using pre-simulated Wiener increments. The user can select
#' different diffusion types via the `model` argument.
#' The update rule is based on Algorithm 8.5 in **Simo Särkkä** and **Arno Solin**: *Applied Stochastic Differential Equations*.

#'
#' @examples
#' # Simulate an Ornstein-Uhlenbeck process
#' simulate_pearson_diffusion(
#'   step_length = 0.1,
#'   par = c(beta = 0.5, mu = -0.5, sigma = 0.1),
#'   total_time = 10,
#'   model = "OU"
#' )
#'
#' # Simulate a square-root process
#' simulate_pearson_diffusion(
#'   step_length = 0.1,
#'   par = c(beta = 0.5, mu = 2, sigma = 0.1),
#'   total_time = 10,
#'   model = "sqrt"
#' )
#'
#' @export
simulate_pearson_diffusion <- function(step_length,
                                       par,
                                       X_0 = NA,
                                       total_time,
                                       sample_method = "auto",
                                       model = "OU"){
  beta  <- par[1]
  mu    <- par[2]
  sigma <- par[3]

  # Initial point in fixed point of process if none specified
  if(is.na(X_0)){
  X_0   <- mu
  }
  N     <- as.integer((total_time) / step_length)

  # Initialize the process

  dW    <- fast_rnorm(N, mean = 0, standard_deviation = sqrt(step_length), method = sample_method)
  time  <- step_length * 0:N

  switch(model,
         "OU" = X_t <- updatestep_OU_process(
           X_0 = X_0,
           beta = beta,
           mu = mu,
           sigma = sigma,
           step_length = step_length,
           dW = dW),
         "sqrt" = X_t <- updatestep_sqrt_process(
           X_0 = X_0,
           beta = beta,
           mu = mu,
           sigma = sigma,
           step_length = step_length,
           dW = dW),
         "linear" = X_t <- updatestep_linear_process(
           X_0 = X_0,
           beta = beta,
           mu = mu,
           sigma = sigma,
           step_length = step_length,
           dW = dW),
         "t-dist" = X_t <- updatestep_t_diffusion_process(
           X_0 = X_0,
           beta = beta,
           mu = mu,
           sigma = sigma,
           step_length = step_length,
           dW = dW),
         "F-dist" = X_t <- updatestep_F_diffusion_process(
           X_0 = X_0,
           beta = beta,
           mu = mu,
           sigma = sigma,
           step_length = step_length,
           dW = dW),
         "Jacobi" = X_t <- updatestep_jacobi_diffusion_process(
           X_0 = X_0,
           beta = beta,
           mu = mu,
           sigma = sigma,
           step_length = step_length,
           dW = dW)
  )

  data.frame(t = time,
             X_t)
}


#' Generic Euler-Maruyama scheme based Simulation for SDE-objects
#'
#' This function simulates a general one-dimensional SDE using the Euler-Maruyama method. It will be
#' used automatically if the SDE object does not provide a custom `simulate_fun`.
#'
#' @param step_length Numeric. Time increment for each simulation step.
#' @param par Numeric vector. Model parameters.
#' @param X_0 Numeric. Initial value of the state variable.
#' @param total_time Numeric. Total simulation time.
#' @param drift Function. Drift function drift(x, par).
#' @param diffusion Function. Diffusion function diffusion(x, par).
#' @param sample_method Character. Method for generating random noise ("auto", "stats", "dqrng").
#'
#' @return A data.frame with columns `t` and `X_t`.
#' @keywords internal
simulate_generic <- function(step_length,
                             par,
                             X_0,
                             total_time,
                             drift,
                             diffusion,
                             sample_method = "auto") {

  # Validate types and values
  if (!is.numeric(step_length) || length(step_length) != 1 || step_length <= 0) {
    stop("step_length must be a single positive numeric value")
  }

  if (!is.numeric(par) || length(par) < 1) {
    stop("par must be a numeric vector of length >= 1")
  }

  if (!is.numeric(X_0) || length(X_0) != 1 || is.na(X_0)) {
    stop("X_0 must be a single numeric value and cannot be NA")
  }

  if (!is.numeric(total_time) || length(total_time) != 1 || total_time <= 0) {
    stop("total_time must be a single positive numeric value")
  }

  if (!is.function(drift)) {
    stop("drift must be a function")
  }

  if (!is.function(diffusion)) {
    stop("diffusion must be a function")
  }

  if (!is.character(sample_method) || !(sample_method %in% c("auto", "stats", "dqrng"))) {
    stop('sample_method must be one of "auto", "stats", or "dqrng"')
  }

  N <- as.integer(total_time / step_length)
  t <- seq(0, total_time, length.out = N + 1)
  x <- numeric(N + 1)

  x[1] <- X_0

  dW <- fast_rnorm(N, mean = 0, standard_deviation = sqrt(step_length), method = sample_method)

  for (i in 1:N) {
    drift_val <- drift(x[i], par)
    diffusion_val <- diffusion(x[i], par)
    x[i + 1] <- x[i] + drift_val * step_length + diffusion_val * dW[i]
  }

  data.frame(t = t, X_t = x)
}

