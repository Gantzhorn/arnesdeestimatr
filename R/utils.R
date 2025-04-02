# Internal helper function for fast normal sampling
# Uses dqrng::dqrnorm() if available for speed; otherwise falls back to stats::rnorm().
# Ensures valid method choice and informs users of corrections.

fast_rnorm <- function(number_of_samples, mean = 0, standard_deviation = 1, method = "auto") {
  valid_methods <- c("auto", "stats", "dqrng")

  if (!method %in% valid_methods) {
    message(paste0("Invalid method '", method, "' provided. Defaulting to 'auto'."))
    method <- "auto"
  }

  if (method == "stats" || (method == "auto" && !requireNamespace("dqrng", quietly = TRUE))) {
    return(stats::rnorm(n = number_of_samples, mean = mean, sd = standard_deviation))
  }

  return(dqrng::dqrnorm(n = number_of_samples, mean = mean, sd = standard_deviation))
}
