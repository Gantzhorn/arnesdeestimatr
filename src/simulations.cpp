#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericVector updatestep_additive_model(double X_0,
                                              Rcpp::NumericVector lambda_t,
                                              double A,
                                              double m,
                                              double sigma,
                                              double step_length,
                                              Rcpp::NumericVector dW){


  int N = lambda_t.size();

  Rcpp::NumericVector X_t(N);

  X_t[0] = X_0;

  for(int i = 0; i < (N - 1); i++){
    X_t[i + 1] = X_t[i] - (A * (X_t[i] - m) * (X_t[i] - m) + lambda_t[i]) * step_length + sigma * dW[i] -
      2 * A * (X_t[i] - m) * sigma * (0.5 * dW[i] * step_length) +
      (A * (X_t[i] - m) * (A * (X_t[i] - m) * (X_t[i] - m) + lambda_t[i]) - 0.5 * A * sigma * sigma) * step_length * step_length;
  }

  return X_t;
}

// [[Rcpp::export]]
Rcpp::NumericVector updatestep_sqrt_model(double X_0,
                                              Rcpp::NumericVector lambda_t,
                                              double A,
                                              double m,
                                              double sigma,
                                              double step_length,
                                              Rcpp::NumericVector dW){


  int N = lambda_t.size();

  Rcpp::NumericVector X_t(N);

  X_t[0] = X_0;

  for (int i = 0; i < N - 1; i++) {
    X_t[i + 1] = X_t[i]
    - (A * (X_t[i] - m) * (X_t[i] - m) + lambda_t[i]) * step_length
    + sigma * std::sqrt(X_t[i]) * dW[i]
    + 0.25 * sigma * sigma * (dW[i] * dW[i] - step_length)
    - A * (X_t[i] - m) * sigma * std::sqrt(X_t[i]) * (dW[i] * step_length)
    + ((A * (X_t[i] - m) * (X_t[i] - m) + lambda_t[i]) * (A * (X_t[i] - m))
    - 0.5 * A * sigma * sigma * X_t[i]) * step_length * step_length
    - (1.0 / (4 * std::sqrt(X_t[i]))) *
    ((A * (X_t[i] - m) * (X_t[i] - m) + lambda_t[i]) + sigma * sigma * sigma / 4.0) *
    (dW[i] * step_length);
  }


  return X_t;
}
