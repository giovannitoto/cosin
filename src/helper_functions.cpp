/* -------------------------------------------------------------------------- */

#include <RcppArmadillo.h>
using namespace Rcpp;

/* -------------------------------------------------------------------------- */

arma::vec join_elem(const arma::vec& v1, const int& v2) {
  return arma::join_cols(v1, v2 * arma::ones(1));
}

arma::vec join_elem(const arma::vec& v1, const double& v2) {
  return arma::join_cols(v1, v2 * arma::ones(1));
}

arma::vec join_elem(const int& v1, const arma::vec& v2) {
  return arma::join_cols(v1 * arma::ones(1), v2);
}

arma::vec join_elem(const double& v1, const arma::vec& v2) {
  return arma::join_cols(v1 * arma::ones(1), v2);
}

/* -------------------------------------------------------------------------- */

arma::mat rbinom_vec(const int& len, const int& size, const double& prob) {
  arma::vec v(len);
  for (int i = 0; i < len; i++) {
    v(i) = R::rbinom(size, prob);
  }
  return v;
}

/* -------------------------------------------------------------------------- */

arma::mat rnorm_mat(const int& rows, const int& cols, const double& mean, const double& sd) {
  arma::mat M(rows, cols);
  for (int j = 0; j < cols; j++) {
    for (int i = 0; i < rows; i++) {
      M(i, j) = R::rnorm(mean, sd);
    }
  }
  return M;
}

arma::mat rnorm_vec(const int& len, const double& mean, const double& sd) {
  arma::vec v(len);
  for (int i = 0; i < len; i++) {
    v(i) = R::rnorm(mean, sd);
  }
  return v;
}

/* -------------------------------------------------------------------------- */

arma::mat runif_mat(const int& rows, const int& cols, const double& minVal, const double& maxVal) {
  arma::mat M(rows, cols);
  for (int j = 0; j < cols; j++) {
    for (int i = 0; i < rows; i++) {
      M(i, j) = R::runif(minVal, maxVal);
    }
  }
  return M;
}

arma::mat runif_vec(const int& len, const double& minVal, const double& maxVal) {
  arma::vec v(len);
  for (int i = 0; i < len; i++) {
    v(i) = R::runif(minVal, maxVal);
  }
  return v;
}

/* -------------------------------------------------------------------------- */
