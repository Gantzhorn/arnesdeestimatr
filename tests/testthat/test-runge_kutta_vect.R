test_that("Runge-kutta can solve with scalar boundary condition, h = 0.2 and n = 10", {
  # ODE IN TEST: y' = x + y, y(b) = a
  # SOLUTION: y(x) = (a + b + 1) * exp(x - b) - x - 1

  x_boundary <- 0

  y_boundary <- 0

  step_size <- 0.2

  num_substeps <- 10

  test_value <- runge_kutta_vect(x0 = x_boundary, y0 = y_boundary,
                                 h = step_size, func = function(x,y){x + y}, n = num_substeps)

  expected_value <- (x_boundary + y_boundary + 1) * exp((step_size + x_boundary) - y_boundary) - (step_size + x_boundary) - 1

  expect_equal(test_value, expected_value, tolerance = 1e-7)
})

test_that("Runge-kutta can solve with vector boundary condition, h = 1, and n = 50", {
  # ODE IN TEST: y' = x + y, y(b) = a
  # EXACT KNOWN SOLUTION: y(x) = (a + b + 1) * exp(x - b) - x - 1

  x_boundary <- seq(from = 0, to = 4, length.out = 5)

  y_boundary <- seq(from = 0, to = 2, length.out = 5)

  step_size <- 1

  num_substeps <- 50

  test_values <- runge_kutta_vect(x0 = x_boundary, y0 = y_boundary,
                                 h = step_size, func = function(x,y){x + y}, n = num_substeps)

  expected_values <- (x_boundary + y_boundary + 1) * exp(step_size) - (step_size + x_boundary) - 1

  expect_equal(test_values, expected_values)
})

test_that("Runge-kutta can broadcast x0 onto y0 and vice-versa", {
  # ODE IN TEST: y' = x + y, y(b) = a
  # EXACT KNOWN SOLUTION: y(x) = (a + b + 1) * exp(x - b) - x - 1

  x_scalar <- 0

  x_vector <- seq(from = 0, to = 4, length.out = 5)

  y_scalar <- 1

  y_vector <- seq(from = 0, to = 4, length.out = 5)

  step_size <- 1

  num_substeps <- 1

  x_broadcast <- runge_kutta_vect(x0 = x_scalar, y0 = y_vector,
                                  h = step_size, func = function(x,y){x + y}, n = num_substeps)

  x_full_vector <- runge_kutta_vect(x0 = rep(x_scalar, times = 5), y0 = y_vector,
                                    h = step_size, func = function(x,y){x + y}, n = num_substeps)

  expect_equal(x_broadcast, x_full_vector)

  y_broadcast <- runge_kutta_vect(x0 = x_vector, y0 = y_scalar,
                                  h = step_size, func = function(x,y){x + y}, n = num_substeps)

  y_full_vector <- runge_kutta_vect(x0 = x_vector, y0 = rep(y_scalar, times = 5),
                                    h = step_size, func = function(x,y){x + y}, n = num_substeps)

  expect_equal(x_broadcast, x_full_vector)
})
