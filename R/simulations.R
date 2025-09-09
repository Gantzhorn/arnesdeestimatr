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
#' @param tau Numeric. Time scale over which the bifurcation parameter changes. Default is 100.
#' @param t_0 Numeric. Initial time offset before the bifurcation starts. Default is 10.
#' @param X_0 Numeric. Initial value of the state variable. If `NA`, the function sets it to the stable fixed point. Default is `NA`.
#' @param beyond_tipping Numeric. Additional simulation time beyond the bifurcation point. Default is 0.
#' @param sample_method Character. Method used for generating random noise. Passed to `fast_rnorm()`. Default is `"auto"` other valid arguments are: `"stats"` and `"dqrng"`.
#' @param noise_term Character. Type of stochastic noise to use. Options are:
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
                                                  tau = 100,
                                                  t_0 = 10,
                                                  X_0 = NA,
                                                  beyond_tipping = 0,
                                                  sample_method = "auto",
                                                  noise_term = "additive"){
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

  switch(noise_term,
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

#' @export
simulate_pearson_diffusion <- function(step_length,
                                       par,
                                       total_time,
                                       X_0 = NA,
                                       sample_method = "auto",
                                       model = "OU"){
  beta  <- par[1]
  mu    <- par[2]
  sigma <- par[3]

  # Initial point in fixed point of process if none specified
  if(is.na(X_0)){X_0 <- mu}
  N          <- as.integer((total_time) / step_length)

  # Initialize the process

  dW            <- fast_rnorm(N, mean = 0, standard_deviation = sqrt(step_length), method = sample_method)
  time          <- step_length * 0:N

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
           dW = dW)
  )

  data.frame(t = time,
             X_t)
}

