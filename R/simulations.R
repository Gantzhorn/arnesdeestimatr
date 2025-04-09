simulate_stochastic_saddlenode_bifurcation <- function(step_length, par,
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
         "additive" = X_t <- updatestep_additive_model(
                                X_0 = X_0,
                                lambda_t = lambda_t,
                                A = A,
                                m = m,
                                sigma = sigma,
                                step_length = step_length,
                                dW = dW),
         "sqrt" = X_t <- updatestep_sqrt_model(
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



