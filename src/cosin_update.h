#ifndef cosin_update_H
#define cosin_update_H

arma::mat truncnorm_lg(const arma::mat& y_lower, const arma::mat& y_upper,
                       const arma::mat& mu, const arma::vec& sigma,
                       const arma::mat& u_rand);

arma::mat update_gammaT(const arma::mat& wT, const arma::mat& Beta,
                       const double& prec_beta, const double& prec_gammaT);

arma::mat update_beta(const int& j, const arma::mat& Qbet, const arma::mat& x,
                      const arma::mat& Z_res, const arma::vec& ps, const arma::mat& GammaT,
                      const arma::mat& wT);

arma::mat update_eta(const arma::mat& Lambda, const arma::vec& ps, const arma::mat& Z);

arma::mat update_gammaB(const int& h, const arma::mat& wB, const arma::mat& Dt,
                        const arma::mat& Bh_1, const arma::mat& Phi_L);

arma::mat update_Lambda_star(const int& j, const arma::mat& etarho, const arma::mat& Phi,
                             const arma::mat& Plam, const arma::vec& ps, const arma::mat& Z);

int update_d(const int& h, const arma::mat& Phi, const arma::vec& rho,
             const arma::mat& eta, const arma::mat& lambdastar,
             const arma::mat& Z, const arma::vec& ps, const arma::vec& w);

arma::mat update_Phi(arma::mat& Phi, const arma::vec& rho, const arma::mat& logit,
                     const double& p_constant, const arma::mat& eta,
                     const arma::mat& lambdastar, const arma::mat& Z, const arma::vec& ps);

#endif
