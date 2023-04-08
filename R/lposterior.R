# ---------------------------------------------------------------------------- #

#' Log-posterior probability
#'
#' Compute the log-posterior probabilities of part or all the MCMC iterations.
#' To speed up computations, it is possible to evaluate only a random fraction of the iterations, using the \code{frac_sampled} argument, or to define the iterations of interest, using the \code{samples} argument.
#'
#' @param out_MCMC A list containing the results of an Adaptive Gibbs Sampler obtained using the function \code{\link{cosin}}.
#' @param frac_sampled A value in (0,1] specifying the fraction of the iterations to evaluate. Default is 1, that is all the iterations.
#' @param samples A vector of integers between 1 and the length of the MCMC chain which specifies the iterations to evaluate. Default is all the iterations.
#' @param columns A string specifying whether to consider only active factors (\code{"kstar"}) or both active and inactive factors (\code{"k"}). Default is \code{"k"}.
#' @param parameters A vector containing the names of the parameters for which you want to compute the posterior mean. The possible valid strings are \code{"beta"}, \code{"gammaT"}, \code{"gammaB"}, \code{"eta"}, \code{"lambda"} and \code{"sigmacol"}. Default is \code{all}, which is equivalent to writing \code{c("beta", "gammaT", "gammaB", "eta", "lambda", "sigmacol")}.
#' @param seed Seed. Default is 28.
#'
#' @return A list with the following elements:
#' \itemize{
#' \item \code{sampled}: indices of the iterations evaluated by the function.
#' \item \code{lposterior}: log-posterior probabilities of the iterations evaluated by the function.
#' \item \code{lposterior_max}: maximum log-posterior probability.
#' \item \code{iteration_max}: index of the iteration in which the log-posterior probability reaches its maximum.
#' \item \code{beta_max}: value of \eqn{\beta} at \code{iteration_max}.
#' \item \code{gammaT_max}: value of \eqn{\Gamma_T} at \code{iteration_max}.
#' \item \code{gammaB_max}: value of \eqn{\Gamma_B} at \code{iteration_max}.
#' \item \code{eta_max}: value of \eqn{\eta} at \code{iteration_max}.
#' \item \code{lambda_max}: value of \eqn{\Lambda} at \code{iteration_max}.
#' \item \code{sigmacol_max}: value of \eqn{\sigma^{-2}_1,\ldots,\sigma^{-2}_p} at \code{iteration_max}.
#' }
#' If \code{columns="k"}, \code{gammaB_max}, \code{eta_max} and \code{lambda_max} contain only columns linked to active factors; if \code{columns="kstar"}, they contain columns linked to both active and inactive factors.
#'
#' @seealso This function is applied to an output of \code{\link{cosin}}.
# The function \code{\link{lposterior_function}} is used to compute the log-posterior probability of a single MCMC iteration.
#'
#' @export
lposterior <- function(out_MCMC, frac_sampled = 1, samples = NULL,
                       columns = "k", parameters = "all", seed = 28) {
  # set seed
  if((length(seed) != 1) || !is.numeric(seed) || (seed != round(seed))) {
    stop("'seed' not valid: it must be an integer.")
  } else {
    set.seed(seed)
  }
  # check whether all necessary variables are available in out_MCMC
  required_variables <- c("y", "wT", "wB", "numFactors", "gammaB", "eta",
                          "lambda", "sigmacol", "hyperparameters")
  if(check_list(out_MCMC, required_variables) == FALSE) {
    stop(paste(paste(c(required_variables), collapse = ", "), "must be stored in the MCMC output."))
  }
  required_hyperparameters <- c("alpha", "a_theta", "b_theta", "sd_gammaT", "sd_beta",
                                "sd_gammaB","a_sigma", "b_sigma", "p_constant", "y_max")
  if(check_list(out_MCMC$hyperparameters, required_hyperparameters) == FALSE) {
    stop(paste(paste(c(required_hyperparameters), collapse = ", "), "must be stored in 'out_MCMC$hyperparameters'."))
  }
  if("x" %in% names(out_MCMC)) {
    if(check_list(out_MCMC, c("beta", "gammaT")) == FALSE) {
      stop("If 'x' is stored in the MCMC output, then also 'beta', 'gammaT' must be stored in the MCMC output.")
    }
    if(check_list(out_MCMC$hyperparameters, c("sd_gammaT", "sd_beta")) == FALSE) {
      stop("If 'x' is stored in the MCMC output, then also 'sd_gammaT', 'sd_beta' must be stored in the hyperparameters of the MCMC output.")
    }
  }
  # number of iterations
  t <- length(out_MCMC$numFactors)
  # sampling
  if((length(frac_sampled) != 1) || (frac_sampled <= 0) || (frac_sampled > 1)) {
    stop("'frac_sampled' not valid: it must be a number in (0,1].")
  }
  if(is.null(samples)) {
    sampled <- sample(1:t, ceiling(t * frac_sampled))
  } else {
    if(any(sapply(samples, function(x) x != round(x)))) {
      stop("'samples' not valid: it must be a vector of integers.")
    }
    if((min(samples) <= 0) || (max(samples) > t)) {
      stop("'samples' not valid: it must be a vector of integers between 1 and the number of iterations.")
    }
    sampled <- sample(samples, ceiling(length(samples) * frac_sampled))
  }
  # sort 'sampled' into ascending order
  sampled <- sort(sampled)
  # check the value of 'columns'
  if((length(columns) != 1) || !(columns %in% c("k", "kstar"))) {
    stop("'columns' not valid: it must be 'k' or 'kstar'.")
  }
  # remove invalid strings from 'parameters'
  valid_parameters_k <- c("gammaB", "eta", "lambda")   # depend on 'columns' argument
  valid_parameters <- c("beta", "gammaT", "sigmacol")  # do not depend on 'columns' argument
  #valid_parameters <- c(valid_par_k, valid_par)
  if("all" %in% parameters) {
    parameters_k <- valid_parameters_k
    parameters <- valid_parameters
  } else {
    parameters_k <- intersect(parameters, valid_parameters_k)
    parameters <- intersect(parameters, valid_parameters)
  }
  # -------------------------------------------------------------------------- #
  # compute the log-posterior probability of each iteration of a MCMC chain
  lpost <- unlist(lapply(sampled, lposterior_function, out_MCMC = out_MCMC, columns = columns))
  # prepare output
  output <- list(sampled = sampled,
                 lposterior = lpost,
                 lposterior_max = max(lpost),
                 iteration_max = sampled[which.max(lpost)]
                )
  # add gammaB, eta, lambda to 'output'
  for (par in parameters_k) {
    if(columns == "k") {
      output[[paste(par, "max", sep="_")]] <- out_MCMC[[par]][[output$iteration_max]]
    } else if(columns == "kstar") {
      output[[paste(par, "max", sep="_")]] <- out_MCMC[[par]][[output$iteration_max]][, 1:out_MCMC$numFactors[output$iteration_max]]
    }
  }
  # add beta, gammaT and sigmacol to 'output'
  for (par in parameters) {
    output[[paste(par, "max", sep="_")]] <- out_MCMC[[par]][[output$iteration_max]]
  }
  return(output)
  # -------------------------------------------------------------------------- #
}

# ---------------------------------------------------------------------------- #

#' Log-posterior probability of a single MCMC iteration
#'
#' Compute the log-posterior probability of a single MCMC iteration.
#'
#' @param ind An integer specifying the index of the iteration of interest.
#' @param out_MCMC A list containing the results of an Adaptive Gibbs Sampler obtained using the function \code{\link{cosin}}.
#' @param columns a string specifying whether to consider only active factors (\code{kstar}) or both active and inactive factors (\code{k}). Default is \code{"k"}.
#'
#' @return A number which is the log-posterior probability of the \code{ind}th iteration of a MCMC chain.
#'
#' @seealso This function is applied to an output of \code{\link{cosin}}. This function is used in \code{\link{lposterior}}.
#'
#' @importFrom stats dgamma dnorm plogis rnorm
lposterior_function <- function(ind, out_MCMC, columns = "k") {
  hyperpar <- out_MCMC$hyperparameters
  # dimension of the data
  n <- nrow(out_MCMC$y)
  p <- ncol(out_MCMC$y)
  # max number of factors (kmax)
  K <- max(out_MCMC$numFactors)
  # complete lambda, beta and eta with the correct number of factors
  rpar <- rescale_parameter(out_MCMC, parameters = "all", columns = columns, samples = ind)
  GammaB <- rpar$gammaB[[ind]]
  eta <- rpar$eta[[ind]]
  Lambda <- rpar$lambda[[ind]]
  rm(rpar)
  # log-likelihood for the model
  if("x" %in% names(out_MCMC)) {
    mean_pr <- tcrossprod(out_MCMC$x, out_MCMC$beta[[ind]])
  } else {
    mean_pr <- matrix(0, n, p)
  }
  # generate values from a truncated normal distribution
  mean_z <- mean_pr + tcrossprod(eta, Lambda)
  Z <- gen_z(Y = out_MCMC$y, mean = mean_z, ps = out_MCMC$sigmacol[[ind]], y_max = hyperpar$y_max)
  # compute inverse and log determinant of omega
  omega_inv <- omega_inversion(lambda = Lambda, sigmai = out_MCMC$sigmacol[[ind]])
  omega_log_det <- det_log_omega(lambda = Lambda, sigmai = out_MCMC$sigmacol[[ind]])
  # log multinorm
  val <- rep(NA, n)
  for(i in 1:n) {
    val[i] <- log_multinorm(z = Z[i, ], mu = mean_pr[i, ],
                            omega_inv = omega_inv,
                            omega_log_det = omega_log_det)
  }
  ll <- sum(val)
  # prior of Lambda beta
  p_gamma1 <- (hyperpar$alpha / (1 + hyperpar$alpha)) ^ seq(1, K)
  prob <- plogis(out_MCMC$wB %*% GammaB)
  p_phi1 <- prob * hyperpar$p_constant  # p_phi1 <-  prob * 2*exp(1)*log(p)/p
  # log prior of beta
  lp_gammaB <- sum(dnorm(GammaB, mean = 0, sd = sqrt(hyperpar$sd_gammaB), log = TRUE))
  lp_col_lambda <- rep(0, K)
  for (h in 1:K) {
    w0  <- which(Lambda[, h] == 0)
    w0l <- length(w0)
    if(w0l == p) {
      lp1 <- log(p_gamma1[h]) + sum(log(1 - p_phi1[,h]))
      lp_col_lambda[h] <- matrixStats::logSumExp(c(log(1 - p_gamma1[h]),  lp1))
    } else if(w0l > 0) {
      lp_col_nonzero <- LaplacesDemon::dmvt(x = Lambda[-w0, h], mu = rep(0, p - w0l),
                                            S = (hyperpar$b_theta / hyperpar$a_theta) * diag(p - w0l),
                                            df = 2 * hyperpar$a_theta, log = TRUE)
      lp_col_lambda[h] <- log(p_gamma1[h]) + sum(log(1 - p_phi1[w0, h])) + sum(log(p_phi1[-w0, h])) + lp_col_nonzero
    } else {
      lp_col_nonzero <- LaplacesDemon::dmvt(x = Lambda[, h], mu = rep(0, p),
                                            S = (hyperpar$b_theta / hyperpar$a_theta) * diag(p),
                                            df = 2 * hyperpar$a_theta, log = TRUE)
      lp_col_lambda[h] <- log(p_gamma1[h]) + sum(log(p_phi1[, h])) + lp_col_nonzero
    }
  }
  lp_lambda <- sum(lp_col_lambda)
  # log prior beta and gammaT
  if("x" %in% names(out_MCMC)){
    # log prior of gammaT
    lp_gammaT = sum(dnorm(out_MCMC$gammaT[[ind]], mean = 0, sd = sqrt(hyperpar$sd_gammaT), log = TRUE))
    # log prior of mu
    me_beta = out_MCMC$wT %*% out_MCMC$gammaT[[ind]]
    lp_beta = sum(dnorm(out_MCMC$beta[[ind]], mean = me_beta, sd = sqrt(hyperpar$sd_beta), log = TRUE))
  } else {
    lp_beta <- lp_gammaT <- 0
  }
  # log prior of sigma (inverse gamma)
  lp_sigma <- sum(dgamma(out_MCMC$sigmacol[[ind]], hyperpar$a_sigma, hyperpar$b_sigma, log = TRUE))
  # log posterior
  lpost <- ll + lp_lambda + lp_gammaB + lp_beta + lp_gammaT + lp_sigma
  return(lpost)
}

# ---------------------------------------------------------------------------- #
# UTILITY FUNCTIONS
# ---------------------------------------------------------------------------- #

#' Logarithm of the determinant of a matrix
#'
#' Compute the logarithm of the determinant of \eqn{\Omega} using the following lemma
#' \deqn{\det(A+UV^\top) = \det(I_m+V^\top A^{-1}U)\det(A)}
#' where \eqn{A} is a \eqn{n\times n} matrix and \eqn{U}, \eqn{V} are \eqn{n\times m} matrices.
#' Since \eqn{\Omega=\Lambda\Lambda^\top+\Sigma}, \eqn{\det(\Omega)} can be written as follows
#' \deqn{\det({\Omega}) = \det(\Lambda\Lambda^\top+\Sigma) = \det(I_k+\Lambda^\top\Sigma^{-1}\Lambda)\det(\Sigma)}
#' where \eqn{\Lambda} is a \eqn{p\times k} matrix and \eqn{\Sigma=\text{diag}(\sigma^2_1,\ldots,\sigma^2_p)} is a \eqn{p\times p} diagonal matrix.
#' Then, \eqn{\log(\det(\Omega))} can be written as follows
#' \deqn{\log(\det(\Omega)) = \log(\det(I_k+\Lambda^\top\Sigma^{-1}\Lambda))+\log(\sigma^2_1)+\ldots+\log(\sigma^2_p)}
#'
#' @param lambda A pxk matrix of factorial weights, \eqn{\Lambda}.
#' @param sigmai A p-dimensional vector containing \eqn{(\sigma^2_1,\ldots,\sigma^2_p)}.
#'
#' @return A number.
#'
#' @seealso This function is used in \code{\link{lposterior_function}}.
det_log_omega <- function(lambda, sigmai) {
  k <- dim(lambda)[2]
  det_log_omega <- log(det(diag(k) + crossprod(lambda, diag(sigmai)) %*% lambda)) + sum(log(sigmai))
  return(det_log_omega)
}

# ---------------------------------------------------------------------------- #

#' Logarithm of a multivariate normal density
#'
#' Compute the logarithm of a multivariate normal density.
#'
#' @param z A vector of values on which density is to be calculated.
#' @param mu A vector of the means of the multivariate normal distribution.
#' @param omega_inv Inverse of the variance matrix of the multivariate normal distribution.
#' @param omega_log_det Logarithm of the determinant of the variance matrix of the multivariate normal distribution.
#'
#' @return A vector with the same dimensions as \code{z}.
#'
#' @seealso This function is used in \code{\link{lposterior_function}}.
log_multinorm <- function(z, mu, omega_inv, omega_log_det) {
  p <- dim(omega_inv)[1]
  zmu <- z - mu
  val <- - p / 2 * log(2 * pi) - 0.5 * omega_log_det - 0.5 * crossprod(zmu, omega_inv) %*% zmu
  return(val)
}

# ---------------------------------------------------------------------------- #

#' Generating values from a truncated normal distribution
#'
#' Generate values from a truncated normal distribution.
#'
#' @param Y A nxp matrix.
#' @param mean A nxp matrix containing the means of the truncated normal distribution.
#' @param ps DA SCRIVERE
#' @param y_max A number specifying the maximum value an element of \code{Y} can take. \code{Inf} is allowed.
#'
#' @return A nxp matrix.
#'
#' @seealso This function is used in \code{\link{lposterior_function}}.
#'
#' @importFrom stats runif
gen_z <- function(Y, mean, ps, y_max) {
  n <- nrow(Y)
  p <- ncol(Y)
  # upper and lower bound
  # rounding function (floor) and the corresponding intervals:
  a_j <- function(j, y_max) {
    val <- j
    val[j == y_max + 1] <- Inf
    val
  }
  # Bounds for truncated normal
  a_y <- a_yp1 <- matrix(NA, nrow = n, ncol = p)
  for (j in 1:p) {
    a_y[,j] <- a_j(Y[, j], y_max)        #   a_y = a_j(y)
    a_yp1[,j] <- a_j(Y[, j] + 1, y_max)  # a_yp1 = a_j(y + 1)
  }
  n_unif <- matrix(runif(n * p), nrow = n, ncol = p)
  Z <- truncnorm_lg(y_lower = log(a_y), y_upper = log(a_yp1),
                    mu = mean, sigma = 1 / sqrt(ps), u_rand = n_unif)
  return(Z)
}

# ---------------------------------------------------------------------------- #
