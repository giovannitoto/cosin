#include <RcppArmadillo.h>
#include "cosin_update.h"      /* truncnorm_lg + update_... */
#include "helper_functions.h"
#include "rcpp_pgdraw.h"       /* samplepg ('pgdraw' R package) */
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// Adaptive Gibbs Sampler (AGS)
// Implementation in C++ of the Adaptive Gibbs Sampler (AGS) for a Generalized Infinite Factor model with Structured Increasing Shrinkage (SIS) prior.
//
// [[Rcpp::export]]
Rcpp::List Rcpp_cosin(double alpha, double a_sigma, arma::mat a_y, arma::mat a_yp1, double a_theta,
          arma::mat GammaB, double b_sigma, double b0, double b1, int burn, arma::mat GammaT, double b_theta,
          arma::vec d,
          arma::mat eta,
          int kmax, int kstar,
          arma::mat Lambda, arma::mat Lambda_star, arma::mat logit,
          arma::mat Beta,
          int nrun,
          Rcpp::List out,
          arma::mat Phi, arma::mat Plam, double prec_gammaT, double prec_beta, arma::mat pred, arma::vec prob, arma::vec ps, double p_constant,
          arma::vec rho,
          double sd_gammaB, int sp, int start_adapt,
          int thin,
          arma::vec uu,
          arma::vec v, bool verbose,
          arma::vec w, arma::mat x, bool xnull,
          arma::mat wT, arma::mat wB) {
  // ---------------------------------------------------------------------------
  // output
  Rcpp::List BETA(sp);
  Rcpp::List GAMMAT(sp);
  Rcpp::List GAMMAB(sp);
  Rcpp::List ETA(sp);
  Rcpp::List LAMBDA(sp);
  Rcpp::List SIG(sp);
  arma::vec K(sp);
  // ---------------------------------------------------------------------------
  // matrix dimensions
  int dd = x.n_cols;
  int k = Lambda_star.n_cols;
  int n = a_y.n_rows;
  int p = a_y.n_cols;
  // int qT = wT.n_cols;
  int qB = wB.n_cols;
  // ---------------------------------------------------------------------------
  int it, i, j, h;
  int ind = 0;
  for (it = 0; it < nrun; it++) {
    if(verbose && (it + 1) % 50 == 0) {
      Rcout << it + 1 << " : " << k << " active factors\n";
    }
    // -------------------------------------------------------------------------
    // 1 - update Z
    arma::mat Zmean(n, p);
    if(xnull) {
      Zmean = eta * trans(Lambda);
    } else {
      Zmean = eta * trans(Lambda) + x * trans(Beta);
    }
    arma::mat n_unif = runif_mat(n, p, 0, 1);
    arma::mat Z = truncnorm_lg(log(a_y), log(a_yp1), Zmean, sqrt(1 / ps), n_unif);
    // -------------------------------------------------------------------------
    if(!xnull) {
      // -----------------------------------------------------------------------
      // 2 - update GammaT
      GammaT = update_gammaT(wT, Beta, prec_beta, prec_gammaT);
      // -----------------------------------------------------------------------
      // 3 - update Beta
      arma::mat Z_res = Z - eta * trans(Lambda);
      arma::mat Qbet(dd, dd);
      if(dd > 1) {
        Qbet = arma::diagmat(prec_beta * arma::ones(dd)) + trans(x) * x;
      } else {
        Qbet = prec_beta + trans(x) * x;
      }
      for (j = 0; j < p; j++) {
        Beta.row(j) = update_beta(j, Qbet, x, Z_res, ps, GammaT, wT);
      }
      Z = Z - x * trans(Beta);
      // -----------------------------------------------------------------------
    }
    // 4 update eta
    eta = update_eta(Lambda, ps, Z);
    // -------------------------------------------------------------------------
    // 5 - update Sigma
    arma::mat Z_res = Z - eta * trans(Lambda);
    for (j = 0; j < p; j++) {
      ps(j) = R::rgamma(a_sigma + 0.5 * n, 1 / (b_sigma + 0.5 * arma::accu(arma::pow(Z_res.col(j), 2))));
    }
    // -------------------------------------------------------------------------
    // 6 - update GammaB
    pred = wB * GammaB;
    logit = arma::exp(pred) / (1 + arma::exp(pred));
    // 6.1 Update phi_L
    arma::mat Phi_L = arma::ones(p, k);
    arma::uvec Phi0 = arma::find(Phi == 0);
    arma::vec logit_phi0 = logit.elem(Phi0);
    arma::uvec which_zero = arma::randu(logit_phi0.n_elem) < (1 - logit_phi0) / (1 - logit_phi0 * p_constant);
    Phi_L.elem(Phi0.elem(arma::find(which_zero))) -= 1;
    // 6.2 Polya gamma
    arma::mat Dt(p, k);
    for (j = 0; j < k; j++) {
      for (i = 0; i < p; i++) {
        Dt(i, j) = samplepg(pred(i, j));
      }
    }
    // 6.3 Update gammaB
    arma::mat Bh_1 = arma::diagmat(arma::ones(qB) / pow(sd_gammaB, 2));
    for (h = 0; h < k; h++) {
      GammaB.col(h) = update_gammaB(h, wB, Dt, Bh_1, Phi_L);
    }
    // -------------------------------------------------------------------------
    // 7 - update Lambda_star and Lambda
    arma::mat etarho = trans(eta.each_row() % trans(rho));
    for (j = 0; j < p; j++) {
      Lambda_star.row(j) = update_Lambda_star(j, etarho, Phi, Plam, ps, Z);
    }
    Lambda = (Lambda_star.each_row() % trans(sqrt(rho))) % Phi;
    // -------------------------------------------------------------------------
    // 8.1 - update d
    for (h = 0; h < k; h++) {
      d(h) = update_d(h, Phi, rho, eta, Lambda_star, Z, ps, w);
    }
    rho = arma::ones(k);
    rho.elem(find(d <= arma::linspace<arma::vec>(0, k - 1, k))) -= 1;
    // 8.2
    arma::vec Plam_diag(k);
    for (h = 0; h < k; h++) {
      Plam_diag(h) = R::rgamma(a_theta + 0.5 * p, 1 / (b_theta + 0.5 * arma::accu(arma::pow(Lambda_star.col(h), 2))));
    }
    Plam = arma::diagmat(Plam_diag);
    // 8.3
    for (h = 0; h < k - 1; h++) {
      v(h) = R::rbeta(1 + arma::accu(d == h), alpha + arma::accu(d > h));
    }
    v(k - 1) = 1;
    w = v % join_elem(1, arma::cumprod(1 - v.head(k - 1)));
    // -------------------------------------------------------------------------
    // 9 - update Phi
    pred = wB * GammaB;
    logit = arma::exp(pred) / (1 + arma::exp(pred));
    Phi = update_Phi(Phi, rho, logit, p_constant, eta, Lambda_star, Z, ps);
    // -------------------------------------------------------------------------
    // save sampled values (after burn-in period)
    if ((it + 1) % thin == 0 && (it + 1) > burn) {
      if (!xnull) {
        if(out.containsElementNamed("beta")) { BETA[ind] = Beta; }
        if(out.containsElementNamed("gammaT")) { GAMMAT[ind] = GammaT; }
      }
      if(out.containsElementNamed("gammaB")) { GAMMAB[ind] = GammaB; }
      if(out.containsElementNamed("eta")) { ETA[ind] = eta; }
      if(out.containsElementNamed("lambda")) { LAMBDA[ind] = Lambda; }
      if(out.containsElementNamed("sigmacol")) { SIG[ind] = ps; }
      if(out.containsElementNamed("numFactors")) { K[ind] = kstar; }
      ind += 1;
    }
    // -------------------------------------------------------------------------
    // Adaptation
    if (uu(it) < prob(it) && (it + 1) > start_adapt) {
      arma::uvec active = find(d > arma::linspace<arma::vec>(0, k - 1, k));
      int kstar_new = active.n_elem;
      kstar = kstar_new;
      if (kstar < k - 1) {
        // set truncation to kstar and subset all variables, keeping only active columns
        k = kstar + 1;
        eta = arma::join_rows(eta.cols(active), rnorm_vec(n, 0, 1));
        double vartheta_k = R::rgamma(a_theta, 1 / b_theta);
        Plam_diag = join_elem(Plam_diag.elem(active), vartheta_k);
        Plam = arma::diagmat(Plam_diag);
        Lambda_star = arma::join_rows(Lambda_star.cols(active), rnorm_vec(p, 0, sqrt(vartheta_k)));
        Phi = arma::join_rows(Phi.cols(active), rbinom_vec(p, 1, p_constant));
        rho = join_elem(rho.elem(active), 1);
        Lambda = arma::join_rows(Lambda.cols(active), Lambda_star.col(k - 1) % Phi.col(k - 1));
        GammaB = arma::join_rows(GammaB.cols(active), rnorm_vec(qB, 0, sqrt(sd_gammaB)));
        w = join_elem(w.elem(active), 1 - sum(w.elem(active)));
        v = join_elem(v.elem(active), 1);
        d = join_elem(d.elem(active), k - 1);
      } else if (k < kmax) {
      // increase truncation by 1 and extend all variables, sampling from the prior/model
      k += 1;
      eta = arma::join_rows(eta, rnorm_vec(n, 0, 1));
      double vartheta_k = R::rgamma(a_theta, 1 / b_theta);
      Plam_diag = join_elem(Plam_diag, vartheta_k);
      Plam = arma::diagmat(Plam_diag);
      Lambda_star = arma::join_rows(Lambda_star, rnorm_vec(p, 0, sqrt(vartheta_k)));
      Phi = arma::join_rows(Phi, rbinom_vec(p, 1, p_constant));
      rho = join_elem(rho, 1);
      Lambda = arma::join_rows(Lambda, Lambda_star.col(k - 1) % Phi.col(k - 1));
      GammaB = arma::join_rows(GammaB, rnorm_vec(qB, 0, sqrt(sd_gammaB)));
      v(k - 2) = R::rbeta(1, alpha);
      v = join_elem(v, 1);
      w = v % join_elem(1, arma::cumprod(1 - v.head(k - 1)));
      d = join_elem(d, k - 1);
      }
    }
    // -------------------------------------------------------------------------
  }
  // ---------------------------------------------------------------------------
  if (!xnull) {
    if(out.containsElementNamed("beta")) { out["beta"] = BETA; }
    if(out.containsElementNamed("gammaT")) { out["gammaT"] = GAMMAT; }
  }
  if(out.containsElementNamed("gammaB")) { out["gammaB"] = GAMMAB; }
  if(out.containsElementNamed("eta")) { out["eta"] = ETA; }
  if(out.containsElementNamed("lambda")) { out["lambda"] = LAMBDA; }
  if(out.containsElementNamed("sigmacol")) { out["sigmacol"] = SIG; }
  if(out.containsElementNamed("numFactors")) { out["numFactors"] = K; }
  // ---------------------------------------------------------------------------
  return out;
  // ---------------------------------------------------------------------------
}
