#' Create a general SDE object
#'
#' @param name Character. Name of the system that the Stochastic Differential Equation models.
#' @param drift Function. Drift function drift(x, par).
#' @param diffusion Function. Diffusion function diffusion(x, par).
#' @param par Numeric vector. Optional ground-truth model parameters used to draw trajectories from the model.
#' @param drift_dx Function. Optional derivative of drift with respect to x_t
#' @param drift_dt Function. Optional derivative of drift with respect to t. This is 0 for all autonomous processes.
#' @param diffusion_dx Function. Optional derivative of diffusion
#' @param lamperti_transform Function. Optional Lamperti transform that omits any parameters. I.e. the anti-derivative of 1 / diffusion(x)
#' @param simulate_fun Function. Optional custom simulation function.
#' @return An object of class \code{SDE}.
#' @export
SDE <- function(name,
                drift,
                diffusion,
                par = NULL,
                drift_dx = NULL,
                drift_dt = NULL,
                diffusion_dx = NULL,
                lamperti_transform = NULL,
                simulate_fun = NULL) {

  if (!is.character(name)) stop("name must be a character string")
  if (!is.function(drift)) stop("drift must be a function")
  if (!is.function(diffusion)) stop("diffusion must be a function")

  obj <- list(
    name = name,
    par = par,
    nparms = length(par),
    drift = drift,
    diffusion = diffusion,
    drift_dx = drift_dx,
    drift_dt = drift_dt,
    diffusion_dx = diffusion_dx,
    lamperti_transform = lamperti_transform,
    simulate_fun = simulate_fun,
    X_t = numeric(),
    t = numeric(),
    simulate_step_length = NA
  )

  class(obj) <- "SDE"
  obj
}


#' Print an SDE
#'
#' Prints a summary of an object of class \code{SDE}.
#'
#' @param x An object of class \code{SDE}.
#' @param ... Additional arguments (currently ignored). Included for S3 consistency.
#' @export
print.SDE <- function(x, ...) {
  cat("SDE object:", x$name, "\n")
  if (is.numeric(x$par)){
    cat("The model is provided with ground truth parameters.\n")
    cat("Number of parameters:", x$nparms, "\n")
    cat("Parameters:", x$par, "\n")
  } else{
    cat("The model is not provided with ground truth parameters.\n")
  }
  if (!is.null(x$simulate_fun)) cat("Has a custom simulation method.\n")
  invisible(x)
}

#' Simulate an SDE
#'
#' Simulates the trajectory of an SDE object.
#'
#' Depending on whether or not a custom simulation method has been given to `simulate_fun`. This method will either:
#' 1. Use the custom simulation function `simulate_fun` of the `SDE`-object and sample according to it, or
#' 2. Use a generic Euler-Maruyama simulation based on the `drift`- and `diffusion`-functions from the `SDE`-object,
#'  in which case the arguments: `step_length`, `total_time`, and `X_0` **must** be provided in `...`.
#'
#' @param object An object of class \code{SDE}.
#' @param step_length Numeric. Time increment for each simulation step.
#' @param ... Additional arguments passed to the simulation function. For the generic simulation,
#'   this must include:
#'   \itemize{
#'     \item \code{total_time}: Numeric. Total simulation time.
#'     \item \code{X_0}: Numeric. Initial value of the state variable.
#'   }
#'  Additionally, \code{sample_method} can be specified for the generic simulation method.
#'  Although this is optional. The valid parameters for this argument is ("auto", "stats", "dqrng"). Default is "auto".
#' @export
simulate_SDE <- function(object, step_length, ...) {
  args <- list(...)

  if (!inherits(object, "SDE")) stop("object must be of class 'SDE'.")
  if (!is.numeric(object$par)) stop("object must have numeric parameters to simulate.")
  if (!is.null(object$simulate_fun)) {
    sim <- object$simulate_fun(par = object$par, ...)
  } else {
    # Ensure required arguments for generic simulation are provided
    required <- c("total_time", "X_0")
    missing_args <- required[!required %in% names(args)]
    if (length(missing_args) > 0) {
      stop("Missing required arguments for generic simulation: ", paste(missing_args, collapse = ", "))
    }

    sim <- simulate_generic(
      step_length   = step_length,
      total_time    = args$total_time,
      X_0           = args$X_0,
      par           = object$par,
      drift         = object$drift,
      diffusion     = object$diffusion,
      sample_method = if ("sample_method" %in% names(args)) args$sample_method else "auto"
    )
  }
  # Store the simulated trajectory in the object
  object$X_t <- sim$X_t
  object$t <- sim$t
  object$simulate_step_length <- step_length
  return(object)
}




