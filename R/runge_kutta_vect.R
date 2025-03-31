#' Vectorized Runge-Kutta Method (RK4)
#'
#' This function implements the **Runge-Kutta 4th Order (RK4)** method for solving
#' ordinary differential equations (ODEs). It is designed to handle **vectorized**
#' input for initial conditions `x0` and `y0`. The function gives the approximation
#' to the solution of the initial value problem at `x0 + h` for all starting positions `x0`
#' based on the provided function, `func`, a step-size `h`, and number of sub-steps `n`.
#'
#' @param x0 A numeric vector of initial values for the independent variable.
#' @param y0 A numeric vector of initial values for the dependent variable.
#' @param h A numeric scalar that defines the step size for the Runge-Kutta method.
#' @param func A function that describes the derivative in the point some point `(x,y)`. It should take `x0` and `y0` as arguments and be vectorized
#' @param n An integer representing the number of sub-steps to iterate through. Setting n = 1 is typically sufficient in estimation applications.
#'
#' @return A numeric vector of estimated values for the function described `func`, i.e. the ODE, evaluated at the points `x0 + h`.
#' The output will have the same length as `x0` and `y0`.
#'
#' @export
#'
#' @examples
#' # Example of using runge_kutta_vect
#' func <- function(x, y) { x + y }  # Simple ODE: y'(x) = x + y
#' x0 <- c(0, 1)  # Initial condition for x
#' y0 <- c(1, 2)  # Initial condition for y
#' h <- 0.1  # Step size
#' n <- 1  # Number of sub-steps
#' result <- runge_kutta_vect(x0, y0, h, func, n)
#' print(result)
runge_kutta_vect <- function(x0, y0, h, func, n) {

  h <- h / n

  for(i in 1:n) {
    k1 <- func(x0, y0)
    k2 <- func(x0 + 0.5 * h, y0 + 0.5 * h * k1)
    k3 <- func(x0 + 0.5 * h, y0 + 0.5 * h * k2)
    k4 <- func(x0 + h, y0 + h * k3)

    y0 <- y0 + h * (k1 + 2 * k2 + 2 * k3 + k4) / 6
    x0 <- x0 + h
  }
  y0
}
