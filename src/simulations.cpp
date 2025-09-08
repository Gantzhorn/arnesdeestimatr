#include <Rcpp.h>

// [[Rcpp::export]]
Rcpp::NumericVector updatestep_saddlenode_additive_model(double X_0,
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
    X_t[i + 1] = X_t[i] -
        (A * (X_t[i] - m) * (X_t[i] - m) + lambda_t[i]) * step_length
        + sigma * dW[i]

      - 2 * A * (X_t[i] - m) * sigma * (0.5 * dW[i] * step_length)

      + (A * (X_t[i] - m) * (A * (X_t[i] - m) * (X_t[i] - m) + lambda_t[i])
          - 0.5 * A * sigma * sigma) * step_length * step_length;
  }

  return X_t;
}

// [[Rcpp::export]]
Rcpp::NumericVector updatestep_saddlenode_sqrt_model(double X_0,
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

    -(1.0 / (4 * std::sqrt(X_t[i]))) *
        ((A * (X_t[i] - m) * (X_t[i] - m) + lambda_t[i]) + sigma * sigma * sigma / 4.0) *
        (dW[i] * step_length);
  }


  return X_t;
}

// [[Rcpp::export]]
Rcpp::NumericVector updatestep_saddlenode_linear_model(double X_0,
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
    + sigma * X_t[i] * dW[i]

    + 0.5 * sigma * sigma * X_t[i] * (dW[i] * dW[i] - step_length)
        - A * (X_t[i] - m) * sigma * X_t[i] * dW[i] * step_length

    + (A * (X_t[i] - m) * ((X_t[i] - m) * (X_t[i] - m) + lambda_t[i])
        - 0.5 * A * sigma * sigma * X_t[i] * X_t[i]) * step_length * step_length

    - 0.5 * sigma * (A * (X_t[i] - m) * (X_t[i] - m) + lambda_t[i]) * dW[i] * step_length;
  }

  return X_t;
}

// [[Rcpp::export]]
Rcpp::NumericVector updatestep_saddlenode_t_dist(double X_0,
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
      +  sigma * std::sqrt((X_t[i] * X_t[i] + 1)) * dW[i]

    + 0.5 * sigma * sigma * X_t[i] * (dW[i] * dW[i] - step_length)
      - A * (X_t[i] - m) * sigma * std::sqrt((X_t[i] * X_t[i] + 1)) * dW[i] * step_length

    + ((A * (X_t[i] - m) * (X_t[i] - m) + lambda_t[i]) * (A * (X_t[i] - m))
      - 0.5 * A * sigma * sigma * (X_t[i] * X_t[i] + 1)) * step_length * step_length

    - 0.5 * (sigma * X_t[i] * (A * (X_t[i] - m) * (X_t[i] - m) + lambda_t[i]) - 0.5 * sigma * sigma * sigma)
        * dW[i] * step_length / std::sqrt((X_t[i] * X_t[i] + 1))
    ;
  }

  return X_t;
}

// [[Rcpp::export]]
Rcpp::NumericVector updatestep_saddlenode_F_dist(double X_0,
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
    +  sigma * std::sqrt(X_t[i] * (X_t[i] + 1)) * dW[i]

    + 0.25 * sigma * sigma * (2 * X_t[i] + 1) * (dW[i] * dW[i] - step_length)
      - A * (X_t[i] - m) * sigma * std::sqrt(X_t[i] * (X_t[i] + 1)) * dW[i] * step_length

    + ((A * (X_t[i] - m) * (X_t[i] - m) + lambda_t[i]) * (A * (X_t[i] - m))
      - 0.5 * A * sigma * sigma * X_t[i] * (X_t[i] + 1)) * step_length * step_length

    - 0.25 * (sigma * (A * (X_t[i] - m) * (X_t[i] - m) + lambda_t[i]) * (2 * X_t[i] + 1) + sigma * sigma * sigma / 4) *
      dW[i] * step_length / sqrt(X_t[i] * (X_t[i] + 1))
    ;
  }

  return X_t;
}


// [[Rcpp::export]]
Rcpp::NumericVector updatestep_saddlenode_jacobi(double X_0,
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
    + sigma * std::sqrt(X_t[i] * (1 - X_t[i])) * dW[i]
    + 0.25 * sigma * (1 - 2 * X_t[i]) * (dW[i] * dW[i] - step_length)

    - A * (X_t[i] - m) * sigma * std::sqrt(X_t[i] * (1 - X_t[i])) * dW[i] * step_length

    + ((A * (X_t[i] - m) * (X_t[i] - m) + lambda_t[i]) * (A * (X_t[i] - m))
         - 0.5 * sigma * sigma * A * (X_t[i] * (1 - X_t[i]))) * step_length * step_length

    - 0.25 * (sigma * (A * (X_t[i] - m) * (X_t[i] - m) + lambda_t[i]) * (1 - 2 * X_t[i]) + sigma * sigma * sigma / 2) *
        dW[i] * step_length / std::sqrt(X_t[i] * (1 - X_t[i]))

    ;
  }

  return X_t;
}
