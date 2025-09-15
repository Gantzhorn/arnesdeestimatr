test_that("SDE constructor validates inputs", {
  drift <- function(x, par) {-par[1] * x}
  diffusion <- function(x, par) {par[2]}

  # valid construction
  sde <- SDE("OU_validates", drift = drift, diffusion = diffusion, par = c(1, 0.5))
  expect_s3_class(sde, "SDE")
  expect_equal(sde$name, "OU_validates")
  expect_equal(sde$nparms, 2)
  expect_equal(sde$par, c(1, 0.5))

  # invalid name
  expect_error(SDE(123, drift, diffusion), "name must be a character string")

  # invalid drift
  expect_error(SDE("bad_drift_fail", drift = 123, diffusion = diffusion), "drift must be a function")

  # invalid diffusion
  expect_error(SDE("bad_diffusion_fail", drift = drift, diffusion = 123), "diffusion must be a function")
})

test_that("simulate_SDE works with fallback Euler-Maruyama", {
  drift <- function(x, par) {-par[1] * x}
  diffusion <- function(x, par) {par[2]}
  sde <- SDE("OU_sim_fallback", drift = drift, diffusion = diffusion, par = c(1, 0.5))

  dqrng::dqset.seed(123)
  sim <- simulate_SDE(sde, step_length = 0.01, total_time = 1, X_0 = 1)

  expect_s3_class(sim, "SDE")
  expect_true(length(sim$X_t) > 2)
  expect_equal(length(sim$t), length(sim$X_t))
  expect_equal(sim$step_length, 0.01)
})

test_that("simulate_SDE works with provided estimation method and gives consistent values to fallback method", {
  drift <- function(x, par) {-par[1] * (x - par[2])}
  diffusion <- function(x, par) {par[3] * sqrt(x)}

  CIR_sim_method <- function(par, step_length, X_0 = NA, total_time, sample_method = "auto"){
    arnesdeestimatr::simulate_pearson_diffusion(step_length, par, X_0, total_time, sample_method, model = "sqrt")
  }

  true_parameters <- c(1, 0.5, 0.1)

  sde_with_method <- SDE("CIR with provided method", drift = drift,
                         diffusion = diffusion,
                         par = true_parameters,
                         simulate_fun = CIR_sim_method)

  sde_without_method <- SDE("CIR with provided method", drift = drift,
                            diffusion = diffusion,
                            par = true_parameters)

  dqrng::dqset.seed(123)
  sde_with_method <- simulate_SDE(sde_with_method, step_length = 0.0001, total_time = 1, X_0 = true_parameters[2])

  dqrng::dqset.seed(123)
  sde_without_method <- simulate_SDE(sde_without_method, step_length = 0.0001, total_time = 1, X_0 = true_parameters[2])

  expect_equal(sde_with_method$X_t, sde_without_method$X_t, tolerance = 0.0001)

})

# TEST OF LIKELIHOOD METHODS

test_that("EM_likelihood_generic returns numeric", {
  drift <- function(x, par) {-par[1] * x}
  diffusion <- function(x, par) {par[2]}
  sde <- SDE("OU_ll_numeric", drift = drift, diffusion = diffusion, par = c(1, 0.5))
  sim <- simulate_SDE(sde, step_length = 0.1, total_time = 1, X_0 = 1)

  ll <- EM_likelihood_generic(sde, par = c(1, 0.5),
                              step_length = sim$step_length,
                              data = sim$X_t)
  expect_true(is.numeric(ll))
  expect_length(ll, 1)
})

test_that("fit_SDE recovers parameters approximately and estimates in expected manner.", {
  drift <- function(x, par) {-par[1] * x}
  diffusion <- function(x, par) {par[2]}
  true_par <- c(5, 1.5)
  sde <- SDE("OU_estimation", drift, diffusion, par = true_par)
  sde <- simulate_SDE(sde, step_length = 0.05, total_time = 100, X_0 = 1)

  fit1 <- fit_SDE(sde, init_par = c(0.5, 1))

  fit2 <- fit_SDE(sde, init_par = c(0.5, 1),
                     step_length = sde$step_length,
                     data = sde$X_t)

  # CHECK THAT THE ONE
  expect_true(is.list(fit1))
  expect_length(fit1$par, 2)
  expect_true(fit1$convergence == 0)


  expect_equal(fit1$par, fit2$par)

})


