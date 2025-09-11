library(testthat)

# --- Tests for simulate_stochastic_saddlenode_bifurcation --- #

test_that("simulate_stochastic_saddlenode_bifurcation gives consistent result for all pearson diffusions", {
  step_length <- 0.1
  par <- c(A = 1.5, m = 0.5, lambda_0 = -0.3, sigma = 0.01)
  t_0 <- 10
  beyond_tipping <- 0
  tau <- 20

  models <- c("additive", "sqrt", "linear", "t-dist", "F-dist", "Jacobi")
  simulationLength <- numeric(length = length(models))

  # ENSURE CONSISTENT STRUCTURE WITHIN EACH MODEL
  for (i in seq_along(models)) {
    model <- models[i]
    res <- simulate_stochastic_saddlenode_bifurcation(
      step_length = step_length,
      par = par,
      t_0 = t_0,
      beyond_tipping = beyond_tipping,
      tau = tau,
      model = model
    )

    simulationLength[i] <- length(res$X_t)

    expect_s3_class(res, "data.frame")
    expect_true(all(c("t", "X_t", "lambda_t", "alpha_t", "mu_t") %in% names(res)))
    expect_equal(length(res$t), length(res$X_t))
    expect_true(is.numeric(res$X_t))
  }

  # ENSURE CONSISTENT STRUCTURE BETWEEN MODEL TYPES
  expect_equal(
    length(unique(simulationLength)), 1
  )

})

test_that("simulate_stochastic_saddlenode_bifurcation throws error for invalid noise term", {
  step_length <- 0.1
  par <- c(A = 1.5, m = 0.5, lambda_0 = -0.3, sigma = 0.01)
  t_0 <- 10
  beyond_tipping <- 0
  tau <- 20

  expect_error(
    simulate_stochastic_saddlenode_bifurcation(
      step_length = step_length,
      par = par,
      t_0 = t_0,
      beyond_tipping = beyond_tipping,
      tau = tau,
      model = "invalid_model"
    )
  )
})


# --- Tests for simulate_pearson_diffusion --- #

test_that("simulate_pearson_diffusion works for all valid models", {
  step_length <- 0.01
  par <- c(beta = 0.75, mu = 0.5, sigma = 0.05)
  total_time <- 10

  models <- c("OU", "sqrt", "linear", "t-dist", "F-dist", "Jacobi")
  simulationLength <- numeric(length = length(models))

  # ENSURE CONSISTENT STRUCTURE WITHIN EACH MODEL
  for (i in seq_along(models)) {
    model <- models[i]
    res <- simulate_pearson_diffusion(
      step_length = step_length,
      par = par,
      total_time = total_time,
      model = model
    )

    simulationLength[i] <- length(res$X_t)

    expect_s3_class(res, "data.frame")
    expect_true(all(c("t", "X_t") %in% names(res)))
    expect_equal(length(res$t), length(res$X_t))
    expect_true(is.numeric(res$X_t))
  }

  # ENSURE CONSISTENT STRUCTURE BETWEEN MODEL TYPES
  expect_equal(
    length(unique(simulationLength)), 1
  )


})

test_that("simulate_pearson_diffusion throws error for invalid model", {
  step_length <- 0.01
  par <- c(beta = 0.75, mu = 0.5, sigma = 0.05)
  total_time <- 10

  expect_error(
    simulate_pearson_diffusion(step_length = step_length,
                               par = par,
                               total_time = total_time,
                               model = "invalid_model")
  )
})
