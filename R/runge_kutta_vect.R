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
