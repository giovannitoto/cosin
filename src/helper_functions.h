#ifndef helper_functions_H
#define helper_functions_H

arma::vec join_elem(const arma::vec& v1, const int& v2);
arma::vec join_elem(const arma::vec& v1, const double& v2);
arma::vec join_elem(const int& v1, const arma::vec& v2);
arma::vec join_elem(const double& v1, const arma::vec& v2);

arma::mat rbinom_vec(const int& len, const int& size, const double& prob);

arma::mat rnorm_mat(const int& rows, const int& cols, const double& mean, const double& sd);
arma::mat rnorm_vec(const int& len, const double& mean, const double& sd);

arma::mat runif_mat(const int& rows, const int& cols, const double& minVal, const double& maxVal);
arma::mat runif_vec(const int& len, const double& minVal, const double& maxVal);

#endif
