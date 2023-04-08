/* -------------------------------------------------------------------------- */

#include <RcppArmadillo.h>
#include <math.h> /* isnan, isinf */
#include "helper_functions.h"
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

/* -------------------------------------------------------------------------- */
// Functions in this file:
//  - truncnorm_lg
//  - update_bmu
//  - update_mu
//  - update_eta
//  - update_beta
//  - update_Lambda_star
//  - update_d
//  - update_Phi
/* -------------------------------------------------------------------------- */

//' Sample from a truncated normal distribution. Samples are drawn
//' componentwise, so each component of the vector is allowed its own
//' mean, standard deviation, and upper and lower limits. The components
//' are assumed to be independent.
//'
//' @param y_lower \code{n x p} matrix of lower endpoints
//' @param y_upper \code{n x p} matrix of upper endpoints
//' @param mu \code{n x p} matrix of conditional expectations
//' @param sigma \code{p x 1} vector of conditional standard deviations
//' @param u_rand \code{n x p} matrix of uniform random variables
//'
//' @return z_star \code{n x p} draw from the truncated normal distribution
//'
//' @note This function uses \code{Rcpp} for computational efficiency.
//'
// [[Rcpp::export]]
arma::mat truncnorm_lg(const arma::mat& y_lower, const arma::mat& y_upper,
                       const arma::mat& mu, const arma::vec& sigma,
                       const arma::mat& u_rand) {
  // Dim of matrix:
  int n = y_lower.n_rows;
  int p = y_lower.n_cols;
  // Storage:
  double val = 0;
  arma::mat z_star(n,p);

  for(int t = 0; t < n; ++t) {
    for(int j = 0; j < p; ++j) {
      // Control
      double uptail1 = (y_lower(t,j) - mu(t,j)) * 1 / sigma(j) > 8;
      double uptail2 = (y_upper(t,j) - mu(t,j)) * 1 / sigma(j) > 8;
      // pnorm(q, mean = 0, sd = 1, lower.tail = TRUE, log.p = FALSE)
      // true + false = 1, true + true = 2, false + false = 0
      if((uptail1 + uptail2) == 0){
        // Lower and upper limits, transformed via pnorm:
        double F_lower = R::pnorm(y_lower(t,j), mu(t,j), sigma(j), 1, 0);
        double F_upper = R::pnorm(y_upper(t,j), mu(t,j), sigma(j), 1, 0);
        // replace 0 with 0.000001 and 1 with 0.999999
        if (F_lower == 0) {
          F_lower = 0.000001;
        } else if (F_upper == 1) {
          F_lower = 0.999999;
        }
        if (F_upper == 0) {
          F_upper = 0.000001;
        } else if (F_upper == 1) {
          F_upper = 0.999999;
        }
        // Corresponding sampled value:
        val = R::qnorm(F_lower + u_rand(t,j) * (F_upper - F_lower), mu(t,j), sigma(j), 1, 0);
      }
      else {
        double F_lower = R::pnorm(y_lower(t,j), mu(t,j), sigma(j), 0, 0);
        double F_upper = R::pnorm(y_upper(t,j), mu(t,j), sigma(j), 0, 0);
        // replace 0 with 0.000001 and 1 with 0.999999
        if (F_lower == 0) {
          F_lower = 0.000001;
        } else if (F_upper == 1) {
          F_lower = 0.999999;
        }
        if (F_upper == 0) {
          F_upper = 0.000001;
        } else if (F_upper == 1) {
          F_upper = 0.999999;
        }
        // Corresponding sampled value:
        val = R::qnorm(F_lower + u_rand(t,j) * (F_upper - F_lower), mu(t,j), sigma(j), 0, 0);
      }
      z_star(t,j) = std::min(std::max(y_lower(t,j), val), y_upper(t,j));
    }
  }
  return z_star;
}

/* -------------------------------------------------------------------------- */

// Update GammaT in the Adaptive Gibbs Sampler
//
// @param wT A pxqT matrix.
// @param Beta A pxd matrix.
// @param prec_beta A positive number.
// @param prec_gammaT A positive number.
//
// @return A qTxd matrix.
//
// @note This function uses \code{Rcpp} for computational efficiency.
arma::mat update_gammaT(const arma::mat& wT, const arma::mat& Beta,
                        const double& prec_beta, const double& prec_gammaT) {
  int d = Beta.n_cols;
  int qT = wT.n_cols;
  arma::mat Trsg = wT * prec_beta;
  arma::mat Vgmbeta1(qT, qT);
  if (qT > 1) {
    Vgmbeta1 = arma::diagmat(prec_gammaT * arma::ones(qT)) + trans(Trsg) * wT;
  } else {
    Vgmbeta1 = prec_gammaT + trans(Trsg) * wT;
  }
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, Vgmbeta1);
  eigval = arma::reverse(eigval);
  eigvec = arma::fliplr(eigvec);
  arma::mat Tmat(qT, qT);
  if (arma::min(eigval) > 0.000001) {
    Tmat = trans(eigvec.each_row() % trans(sqrt(eigval)));
  } else {
    Tmat = arma::chol(Vgmbeta1);
  }
  arma::mat Q;
  arma::mat R;
  arma::qr(Q, R, Tmat);
  arma::mat S = arma::inv(R);
  arma::mat GammaT = trans(trans(Beta) * Trsg * S * trans(S) + rnorm_mat(d, qT, 0, 1) * trans(S));
  return GammaT;
}

/* -------------------------------------------------------------------------- */

// Update the jth row of Beta in the Adaptive Gibbs Sampler
//
// @param j An integer number.
// @param Qbet A dxd matrix.
// @param prec_beta A positive number.
// @param x A nxd matrix.
// @param Z_res A nxp matrix.
// @param ps A p-dimensional vector.
// @param GammaT A dxqT matrix.
// @param wT A pxqT matrix.
//
// @return A 1xd matrix.
//
// @note This function uses \code{Rcpp} for computational efficiency.
arma::mat update_beta(const int& j, const arma::mat& Qbet, const arma::mat& x,
                      const arma::mat& Z_res, const arma::vec& ps, const arma::mat& GammaT,
                      const arma::mat& wT) {
  int d = x.n_cols;
  arma::mat Lbet = trans(arma::chol(Qbet));
  arma::vec bbet = trans(x) * Z_res.col(j) + ps(j) * trans(GammaT) * trans(wT.row(j));
  // mean
  arma::mat vbet = solve(arma::trimatl(Lbet), bbet);
  arma::mat mbet = solve(arma::trimatu(trans(Lbet)), vbet);
  // var
  arma::vec zbet = rnorm_vec(d, 0, 1);
  arma::mat ybet = solve(trimatu(trans(Lbet)), zbet);
  arma::mat betaj = trans(ybet + mbet);
  return betaj;
}

/* -------------------------------------------------------------------------- */

// Update eta in the Adaptive Gibbs Sampler
//
// @param Lambda A pxk matrix.
// @param ps A p-dimensional vector.
// @param Z A nxp matrix.
//
// @return A nxk matrix.
//
// @note This function uses \code{Rcpp} for computational efficiency.
arma::mat update_eta(const arma::mat& Lambda, const arma::vec& ps, const arma::mat& Z) {
  int k = Lambda.n_cols;
  int n = Z.n_rows;
  arma::mat Lmsg = Lambda.each_col() % ps;
  arma::mat Veta1(k, k);
  if (k > 1) {
    Veta1 = arma::diagmat(arma::ones(k)) + trans(Lmsg) * Lambda;
  } else {
    Veta1 = 1 + trans(Lmsg) * Lambda;
  }
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, Veta1);
  eigval = arma::reverse(eigval);
  eigvec = arma::fliplr(eigvec);
  arma::mat Tmat(k, k);
  if (arma::min(eigval) > 0.000001) {
    Tmat = trans(eigvec.each_row() % trans(sqrt(eigval)));
  } else {
    Tmat = arma::chol(Veta1);
  }
  arma::mat Q;
  arma::mat R;
  arma::qr(Q, R, Tmat);
  arma::mat S = arma::inv(R);
  arma::mat eta = Z * Lmsg * S * trans(S) + rnorm_mat(n, k, 0, 1) * trans(S);
  return eta;
}

/* -------------------------------------------------------------------------- */

// Update the hth column of GammaB in the Adaptive Gibbs Sampler
//
// @param h An integer number.
// @param wB A pxqB matrix.
// @param Dt A pxk matrix.
// @param Bh_1 A qBxqB matrix.
// @param Phi_L A pxk matrix.
//
// @return A qBx1 matrix.
//
// @note This function uses \code{Rcpp} for computational efficiency.
arma::mat update_gammaB(const int& h, const arma::mat& wB, const arma::mat& Dt,
                        const arma::mat& Bh_1, const arma::mat& Phi_L) {
  int qB = wB.n_cols;
  arma::mat Qbeta = trans(wB.each_col() % Dt.col(h)) * wB + Bh_1;
  arma::mat Lbeta = trans(arma::chol(Qbeta));
  arma::vec bbeta = trans(wB) * (Phi_L.col(h) - 0.5);
  // mean
  arma::mat vbeta = solve(trimatl(Lbeta), bbeta);
  arma::mat mbeta = solve(trimatu(trans(Lbeta)), vbeta);
  // var
  arma::vec zbeta = rnorm_vec(qB, 0, 1);
  arma::mat ybeta = solve(trimatu(trans(Lbeta)), zbeta);
  arma::mat GammaBh = ybeta + mbeta;
  return GammaBh;
}

/* -------------------------------------------------------------------------- */

// Update the jth row of Lambda_star in the Adaptive Gibbs Sampler
//
// @param j An integer number.
// @param etarho A kxn matrix.
// @param Phi A pxk matrix.
// @param Plam A kxk matrix.
// @param ps A p-dimensional vector;
// @param Z A nxp matrix.
//
// @return A 1xk matrix.
//
// @note This function uses \code{Rcpp} for computational efficiency.
arma::mat update_Lambda_star(const int& j, const arma::mat& etarho, const arma::mat& Phi,
                             const arma::mat& Plam, const arma::vec& ps, const arma::mat& Z) {
  int k = Phi.n_cols;
  arma::mat etaj = trans(etarho.each_col() % trans(Phi.row(j)));
  arma::mat Qlam = Plam + ps(j) * trans(etaj) * etaj;
  arma::mat Llam = trans(arma::chol(Qlam));
  arma::vec blam = ps(j) * trans(etaj) * Z.col(j);
  // mean
  arma::mat vlam = solve(trimatl(Llam), blam);
  arma::mat mlam = solve(trimatu(trans(Llam)), vlam);
  // var
  arma::vec zlam = rnorm_vec(k, 0, 1);
  arma::mat ylam = solve(trimatu(trans(Llam)), zlam);
  arma::mat lambda_starj = trans(ylam + mlam);
  return lambda_starj;
}

/* -------------------------------------------------------------------------- */

// Update the hth element of d in the Adaptive Gibbs Sampler
//
// @param h An integer.
// @param Phi A pxk matrix.
// @param rho A k-dimensional vector.
// @param eta A nxk matrix.
// @param lambdastar A pxk matrix.
// @param Z A nxp matrix.
// @param ps A p-dimensional vector.
// @param w A k-dimensional vector.
//
// @return An integer in 1, ..., k.
//
// @note This function uses \code{Rcpp} for computational efficiency.
int update_d(const int& h, const arma::mat& Phi, const arma::vec& rho,
             const arma::mat& eta, const arma::mat& lambdastar,
             const arma::mat& Z, const arma::vec& ps, const arma::vec& w) {
  int k = Phi.n_cols;
  int n = eta.n_rows;
  int p = Phi.n_rows;
  int i, j, l;
  double lnorm0 = 0.0, lnorm1 = 0.0;
  arma::vec sdy = sqrt(1 / ps);
  for (i = 0; i < n; i++) {
    for (j = 0; j < p; j++) {
      // initialize mean
      double muijh = -1 * sqrt(rho(h)) * Phi(j, h) * lambdastar(j, h) * eta(i, h);
      for (l = 0; l < k; l++) {
        muijh += sqrt(rho(l)) * Phi(j, l) * lambdastar(j, l) * eta(i, l);
      }
      lnorm0 += R::dnorm(Z(i, j), muijh, sdy(j), 1);
      // update mean
      muijh += Phi(j, h) * lambdastar(j, h) * eta(i, h);
      lnorm1 += R::dnorm(Z(i, j), muijh, sdy(j), 1);
    }
  }
  // adjust the scale
  double mlnorm = std::max(lnorm0, lnorm1);
  lnorm0 -= mlnorm;
  lnorm1 -= mlnorm;
  Rcpp::NumericVector prob_h = as<NumericVector>(wrap(log(w)));
  for (i = 0; i <= h; i++) {
    prob_h(i) += lnorm0;
  }
  for (i = h + 1; i < k; i++) {
    prob_h(i) += lnorm1;
  }
  prob_h = exp(prob_h);
  if (sum(prob_h) == 0.0) {
    prob_h = rep(0, k);
    prob_h(k - 1) = 1.0;
  }
  prob_h = prob_h / sum(prob_h);
  // draw d from a multinomial distribution
  Rcpp::IntegerVector d(k);
  R::rmultinom(1, prob_h.begin(), k, d.begin());
  return Rcpp::which_max(d);
}

/* -------------------------------------------------------------------------- */

// Update Phi in the Adaptive Gibbs Sampler
//
// @param rho A k-dimensional vector.
// @param logit A pxk matrix.
// @param p_constant A number in (0,1).
// @param p An integer.
// @param n An integer.
// @param eta A nxk matrix.
// @param lambdastar A pxk matrix.
// @param Phi A pxk matrix.
// @param Z A nxp matrix.
// @param ps A p-dimensional vector.
// @param k An integer.
//
// @return A pxk matrix.
//
// @note This function uses \code{Rcpp} for computational efficiency.
arma::mat update_Phi(arma::mat& Phi, const arma::vec& rho, const arma::mat& logit,
                     const double& p_constant, const arma::mat& eta,
                     const arma::mat& lambdastar, const arma::mat& Z, const arma::vec& ps) {
  // define variables
  int k = Phi.n_cols;
  int n = eta.n_rows;
  int p = Phi.n_rows;
  int i, j, l, h;
  arma::uword f;
  arma::vec p_phi0(p), p_phi1(p), p_phi_sum(p), lnorm0(p), lnorm1(p);
  arma::uvec wr0 = find(rho == 0);
  arma::uvec wr1 = find(rho == 1);
  arma::vec sdy = sqrt(1 / ps);
  // update for inactive factors
  for (f = 0; f < wr0.n_elem; f++) {
    h = wr0(f);
    p_phi0 = 1 - logit.col(h) * p_constant;
    p_phi1 = logit.col(h) * p_constant;
    p_phi_sum = p_phi0 + p_phi1;
    for (j = 0; j < p; j++) {
      Phi(j, h) = R::runif(0, 1) * p_phi_sum(j) < p_phi1(j);
    }
  }
  // update for active factors
  for (f = 0; f < wr1.n_elem; f++) {
    h = wr1(f);
    lnorm0 = arma::zeros(p);
    lnorm1 = arma::zeros(p);
    for (j = 0; j < p; j++) {
      for (i = 0; i < n; i++) {
        // initialize mean
        double muijh = -1 * rho(h) * Phi(j, h) * lambdastar(j, h) * eta(i, h);
        for (l = 0; l < k; l++) {
            muijh += rho(l) * Phi(j, l) * lambdastar(j, l) * eta(i, l);
        }
        lnorm0(j) += R::dnorm(Z(i, j), muijh, sdy(j), 1);
        // update mean
        muijh += rho(h) * lambdastar(j, h) * eta(i, h);
        lnorm1(j) += R::dnorm(Z(i, j), muijh, sdy(j), 1);
      }
      // adjust the scale
      double mlnorm = std::max(lnorm0(j), lnorm1(j));
      lnorm0(j) -= mlnorm;
      lnorm1(j) -= mlnorm;
    }
    p_phi0 = exp(lnorm0 + log(1 - logit.col(h) * p_constant));
    p_phi1 = exp(lnorm1 + log(logit.col(h) * p_constant));
    p_phi_sum = p_phi0 + p_phi1;
    for (j = 0; j < p; j++) {
      Phi(j, h) = R::runif(0, 1) * p_phi_sum(j) < p_phi1(j);
    }
  }
  return Phi;
}

/* -------------------------------------------------------------------------- */
