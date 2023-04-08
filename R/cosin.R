#' AGS for GIF model with COSIN prior
#'
#' @description
#' Implementation of the Adaptive Gibbs Sampler (AGS) for a Generalized Infinite Factor (GIF) model with COunt data Structured INfinite factorization (COSIN) prior.
#'
#' @details
#' Suppose we have a \eqn{n\times p} matrix \eqn{y} of gene expression containing the counts for \eqn{p}
#' genes measured in on \eqn{n} cells, a \eqn{n\times d} matrix \eqn{x} of cell-specific covariates,
#' a \eqn{n\times q_T} matrix \eqn{w_T} of gene-specific technical meta-covariates and
#' a \eqn{n\times q_B} matrix \eqn{w_B} of gene-specific biological meta-covariates.
#'
#' We introduce a continuous-valued latent matriz \eqn{z} related to the observed count-valued matrix \eqn{y}
#' via a STAR operator \eqn{\mathcal{S}:\mathbb{R}\to\mathbb{N}} with \eqn{\mathcal{S}(\cdot)=\mathcal{H}(\mathcal{G}(\cdot))}.
#' In this case, the transformation operator \eqn{\mathcal{G}} is the exponential trasformation and the rounding operator
#' \eqn{\mathcal{H}} is such that \eqn{\mathcal{H}(t)=l} if \eqn{t\in\mathbb{A}_l} and \eqn{\{\mathbb{A}_l\}^\infty_{l=1}}
#' is a known partition of \eqn{\mathbb{R}}.
#' Thus, the single entry \eqn{y_{ij}} of \eqn{y} is linked to a latent \eqn{z_{ij}} of \eqn{z} as
#' \eqn{y_{ij} = l \quad\text{if}\quad\exp\{z_{ij}\}\in[l,l+1)}.
#'
#' The latent variables \eqn{z_{ij}} are modeled via
#' \deqn{z_{ij} = x_i^\top\beta_j+\epsilon_{ij},}
#' \deqn{\epsilon_{ij}=\sum_{h=1}^k\lambda_{jh}\eta_{ih}+\varepsilon_{ij},}
#' \deqn{\varepsilon_{i}\sim N(0,\Sigma),}
#' where \eqn{\beta_j} is a row of the \eqn{p\times} d matrix \eqn{\beta}, \eqn{\lambda_{jh}} is an element of the \eqn{p\times k} factor loading matrix \eqn{\lambda},
#' \eqn{\eta_{ih}} is an element of the \eqn{n\times k} latent factor matrix \eqn{\eta}, and
#' \eqn{\varepsilon_{1},\ldots,\varepsilon_{n}} are iid random errors with \eqn{\Sigma=diag(\sigma^2_1,\ldots,\sigma^2_p)}.
#'
#' Following a Bayesian approach we elicit suitable prior distributions for the model parameters:
#' \deqn{\beta_j\sim N_d(\Gamma_Tw_{jT},\sigma^2_\beta I_d), \quad \gamma_{mT}\sim N_d(0,\sigma^2_{\Gamma_T}I_d),}
#' \deqn{\eta_i\sim N_k(0, I_k), \quad \sigma^{-2}_j\sim Ga(a_\sigma,b_\sigma),}
#' \deqn{\lambda_{jh}\sim N(0,\vartheta_{jh}\rho_h\phi_{jh}), \quad \vartheta^{-1}_h\sim Ga(a_\theta,b_\theta),}
#' \deqn{\rho_h\sim Bern(1-\pi_h), \quad \pi_h=\sum_{l=1}^hv_l\prod_{m=1}^{l-1}v_m, \quad v_m\sim Beta(1,\alpha),}
#' \deqn{E[\phi_{jh}] = c_p\text{logit}^{-1}(w_{jB}^\top\gamma_{hB}), \quad \gamma_{hB}\sim N_{q_B}(0,\sigma^2_{\Gamma_B}I_{q_B}),}
#' where \eqn{\gamma_{mT}} is a column of the \eqn{d\times q_T} coefficient matrix \eqn{\Gamma_T},
#' \eqn{\gamma_{hB}} is a column of the \eqn{q_B\times k} coefficient matrix \eqn{\Gamma_B} and \eqn{c_p\in(0,1)} is a offset.
#'
#' @param y A \eqn{n\times p} matrix \eqn{y} of gene expression.
#' @param wT A matrix \eqn{w_T} of gene-specific technical meta-covariates having \eqn{p} rows; the variables must be numeric or factors.
#' @param wB A matrix \eqn{w_B} of gene-specific biological meta-covariates having \eqn{p} rows; the variables must be numeric or factors.
#' @param x A matrix \eqn{w} of cell-specific covariates having \eqn{n} rows; the variables must be numeric or factors.
#' @param y_max A fixed and known upper bound for the values in \code{y}. Default is \code{Inf}.
#' @param wTformula Formula specifying the gene-specific technical meta-covariates used in the process mean inserted in the model.
#' @param wBformula Formula specifying the gene-specific biological meta-covariates used in the COSIN prior.
#' @param xFormula Formula specifying the covariates inserted in the model; by default all are considered.
#' @param stdwT Logical: if \code{TRUE}, numeric gene-specific technical meta-covariates are standardized; by default they are standardized.
#' @param stdwB Logical: if \code{TRUE}, numeric gene-specific biological meta-covariates are standardized; by default they are standardized.
#' @param stdx Logical: if \code{TRUE}, numeric cell-specific covariates are standardized; by default they are standardized.
#' @param sd_gammaT Standard deviation for \eqn{\gamma_{mT}} prior distribution. Default is \eqn{\sigma_{\Gamma_T}=1}.
#' @param sd_gammaB Standard deviation for \eqn{\gamma_{hB}} prior distribution. Default is \eqn{\sigma_{\Gamma_B}=1}.
#' @param sd_beta Standard deviation for \eqn{\beta_j} prior distribution. Default is \eqn{\sigma_\mu=1}.
#' @param a_theta,b_theta Shape (\code{a_theta}) and rate (\code{b_theta}) parameters for \eqn{\vartheta_{jh}^{-1}}. Default is \eqn{a_{\theta}=b_{\theta}=1}.
#' @param a_sigma,b_sigma Shape (\code{a_sigma}) and rate (\code{b_sigma}) parameters for \eqn{\sigma_j^{-2}} prior distribution. Default is \eqn{a_{\sigma}=b_{\sigma}=1}.
#' @param alpha Non-negative parameter for \eqn{v_m} prior distribution. Default is \eqn{\alpha=5}.
#' @param p_constant Factor probability constant. Default is \eqn{c_p=10e\log(p)/p}.
#' @param kinit An integer minimun number of latent factors. Default is \code{min(floor(log(p)*kval), p)}.
#' @param kmax Maximum number of latent factors. Default is \code{p+1}.
#' @param kval An integer number used to calculate the default value of \code{kinit}. Default is 6.
#' @param nrun An integer number of iterations. Default is 100.
#' @param burn An integer number of burning iterations (number of iterations to discard). Default is \code{round(nrun/4)}.
#' @param thin An integer thinning value (number of iterations to skip between saving iterations). Default is 1.
#' @param start_adapt An integer number of iterations before adaptation. Default is 50.
#' @param b0,b1 Positive constants for the adaptive probability \eqn{p(t)=\exp(-b_0-b_1t)}. Default is \eqn{b_0=1} and \eqn{b_1=5\times 10^{-4}}.
#' @param seed Seed. Default is 28.
#' @param output A vector containing the names of the parameters for which you want to save the draws from the posterior distribution. The possible valid strings are \code{"beta"}, \code{"gammaT"}, \code{"gammaB"}, \code{"eta"}, \code{"lambda"} and \code{"sigmacol"}. Default is \code{"all"}, which is equivalent to writing \code{c("beta", "gammaT", "gammaB", "eta", "lambda", "sigmacol")}.
#' @param verbose Logical: if \code{TRUE}, print the number of active factors every 50 iterations. Default is \code{TRUE}.
#'
#' @return A list with the following elements:
#' \itemize{
#' \item \code{numFactors}: a vector containing the number of active factors at each saved iteration.
#' \item \code{beta}: a list containing draws from the posterior distribution of \eqn{\beta=(\beta_1,\ldots,\beta_p)^\top\in\mathbb{R}^{p\times d}}.
#' \item \code{gammaT}: a list containing draws from the posterior distribution of \eqn{\Gamma_T=(\gamma_{1T},\ldots,\gamma_{q_TT})\in\mathbb{R}^{d\times q_T}}.
#' \item \code{gammaB}: a list containing draws from the posterior distribution of \eqn{\Gamma_B=(\gamma_{1B},\ldots,\gamma_{kB})\in\mathbb{R}^{q_B\times k}}.
#' \item \code{eta}: a list containing draws from the posterior distribution of \eqn{\eta=(\eta_1,\ldots,\eta_n)^\top\in\mathbb{R}^{n\times k}}.
#' \item \code{lambda}: a list containing draws from the posterior distribution of \eqn{\Lambda\in\mathbb{R}^{p\times k}}.
#' \item \code{sigmacol}: a list containing draws from the posterior distribution of \eqn{(\sigma_1^{-2},\ldots,\sigma_p^{-2})\in\mathbb{R}^{p}}.
#' \item \code{time}: total elapsed time.
#' \item \code{y}: the \eqn{n\times p} matrix \eqn{y} of gene expression provided as input to the function.
#' \item \code{wT}: the \eqn{p\times q_T} matrix of gene-specific technical meta-covariates obtained as a result of variable selection, using \code{wTformula}, standardization of numerical variables, using \code{stdwT}, and conversion of factors into dichotomous variables.
#' \item \code{wB}: the \eqn{p\times q_B} matrix of gene-specific biological meta-covariates obtained as a result of variable selection, using \code{wBformula}, standardization of numerical variables, using \code{stdwB}, and conversion of factors into dichotomous variables.
#' \item \code{x}: the \eqn{n\times d} matrix of cell-specific covariates obtained as a result of variable selection, using \code{xFormula}, standardization of numerical variables, using \code{stdx}, and conversion of factors into dichotomous variables.
#' \item \code{hyperparameters}: a list containing the hyperparameters provided as input to the function.
#' }
#'
#' @importFrom stats formula model.matrix plogis rbeta rbinom rgamma rnorm runif
#'
#' @import RcppArmadillo
#'
#' @export
cosin <- function(y, wT = NULL, wB = NULL, x = NULL, y_max = Inf,
                  wTformula = formula("~ ."), wBformula = formula("~ ."), xFormula = formula("~ ."),
                  stdwT = TRUE, stdwB = TRUE, stdx = TRUE,
                  sd_gammaT = 1, sd_gammaB = 1, sd_beta = 1,
                  a_theta = 1, b_theta = 1, a_sigma = 1, b_sigma = 1,
                  alpha = 5, p_constant = NULL,
                  kinit = NULL, kmax = NULL, kval = 6,
                  nrun = 100, burn = round(nrun/4), thin = 1,
                  start_adapt = 50, b0 = 1, b1 = 5*10^(-4),
                  seed = 28, output = "all", verbose = TRUE) {
  # -------------------------------------------------------------------------- #
  # set seed
  if((length(seed) != 1) || !is.numeric(seed) || (seed != round(seed))) {
    stop("'seed' not valid: it must be an integer.")
  } else {
    set.seed(seed)
  }
  # -------------------------------------------------------------------------- #
  if((length(verbose) != 1) || !is.logical(verbose)) {
    stop("'verbose' not valid: it must be 'TRUE' or 'FALSE'.")
  }
  # -------------------------------------------------------------------------- #
  if(!is.data.frame(y) && !is.matrix(y)) {
    stop("'y' not valid: it must be a matrix or a dataframe.")
  }
  # -------------------------------------------------------------------------- #
  n <- dim(y)[1]
  p <- dim(y)[2]
  # -------------------------------------------------------------------------- #
  if(is.null(wT)) {
    wT <- matrix(1, nrow = p, ncol = 1)
  } else {
    if(!is.data.frame(wT) && !is.matrix(wT)) {
      stop("'wT' not valid: it must be a matrix or a dataframe.")
    }
    if(p != nrow(wT)) {
      stop("'y' and 'wT' not compatible: the number of columns of 'y' must be equal to the number of rows of 'wT'.")
    }
    if((length(stdwT) != 1) || !is.logical(stdwT)) {
      stop("'stdwT' not valid: it must be 'TRUE' or 'FALSE'.")
    }
    if(stdwT & is.data.frame(wT)) {
      is.fact.wT = sapply(wT, is.factor)
      wT[, is.fact.wT == FALSE] <- scale(wT[, is.fact.wT == FALSE])
      wT <- model.matrix(wTformula, wT)
    }
  }
  if(is.null(wB)) {
    wB <- matrix(1, nrow = p, ncol = 1)
  } else {
    if(!is.data.frame(wB) && !is.matrix(wB)) {
      stop("'wB' not valid: it must be a matrix or a dataframe.")
    }
    if(p != nrow(wB)) {
      stop("'y' and 'wB' not compatible: the number of columns of 'y' must be equal to the number of rows of 'wB'.")
    }
    if((length(stdwB) != 1) || !is.logical(stdwB)) {
      stop("'stdwB' not valid: it must be 'TRUE' or 'FALSE'.")
    }
    if(stdwB & is.data.frame(wB)) {
      is.fact.wB <- sapply(wB, is.factor)
      wB[, is.fact.wB == FALSE] <- scale(wB[, is.fact.wB == FALSE])
      wB <- model.matrix(wBformula, wB)
    }
  }
  qT <- ncol(wT)
  qB <- ncol(wB)
  # -------------------------------------------------------------------------- #
  xnull <- is.null(x)
  if(!xnull) {
    if(!is.data.frame(x) && !is.matrix(x)) {
      stop("'x' not valid: it must be a matrix or a dataframe.")
    }
    if(n != nrow(x)) {
      stop("'y' and 'x' not compatible: they must have the same number of rows.")
    }
    if((length(stdx) != 1) || !is.logical(stdx)) {
      stop("'stdx' not valid: it must be 'TRUE' or 'FALSE'.")
    }
    if(stdx & is.data.frame(x)) {
      is.fact.x <- sapply(x, is.factor)
      x[, is.fact.x == FALSE] <-  scale(x[, is.fact.x == FALSE])
      x <- model.matrix(xFormula, x)
    }
    d <- ncol(x)
  } else {
    x <- matrix(1, nrow = 1, ncol = 1)  # not used but necessary for Rcpp
    d <- 1
  }
  # -------------------------------------------------------------------------- #
  if((length(kval) != 1) || !is.numeric(kval) || (kval != round(kval))) {
    stop("'kval' not valid: it must be an integer.")
  }
  if(is.null(kinit)) {
    kinit <- min(floor(log(p) * kval), p)
  } else if((length(kinit) != 1) || !is.numeric(kinit) || (kinit != round(kinit))) {
    stop("'kinit' not valid: it must be an integer.")
  }
  if(is.null(kmax)) {
    kmax <- p + 1
  } else if((length(kmax) != 1) || !is.numeric(kmax) || (kmax != round(kmax)) || (kmax < kinit)) {
    stop("'kmax' not valid: it must be an integer greater than or equal to 'kinit'.")
  }
  k <- kinit       # number of factors to start with (active and inactive)
  kstar <- k - 1   # number of active factors
  # -------------------------------------------------------------------------- #
  if((length(nrun) != 1) || !is.numeric(nrun) || (nrun != round(nrun))) {
    stop("'nrun' not valid: it must be an integer.")
  }
  if((length(burn) != 1) || !is.numeric(burn) || (burn != round(burn)) || (burn >= nrun)) {
    stop("'burn' not valid: it must be an integer less than 'nrun'.")
  }
  if((length(thin) != 1) || !is.numeric(thin) || (thin != round(thin)) || (thin > nrun - burn)) {
    stop("'thin' not valid: it must be an integer less than or equal to 'nrun - burn'.")
  }
  if((length(start_adapt) != 1) || !is.numeric(start_adapt) || (start_adapt != round(start_adapt))) {
    stop("'start_adapt' not valid: it must be an integer.")
  }
  # number of posterior samples
  sp <- floor((nrun - burn) / thin)
  # -------------------------------------------------------------------------- #
  if((length(y_max) != 1) || !is.numeric(y_max)) {
    stop("'y_max' not valid: it must be a number or 'Inf'.")
  }
  # Transformation for the response
  a_j <- function(j, y_max) {
    val <- j
    val[j == y_max + 1] <- Inf
    val
  }
  # Bounds for truncated normal
  a_y <- a_yp1 <- matrix(NA, nrow = n, ncol = p)
  for(j in 1:p) {
    a_y[, j] <- a_j(y[, j], y_max)        # a_y = a_j(y)
    a_yp1[, j] <- a_j(y[, j] + 1, y_max)  # a_yp1 = a_j(y + 1)
  }
  # Replace NA with 0/Inf in a_y/a_yp1
  a_y[is.na(y)] <- 0      # log(0)=-Inf
  a_yp1[is.na(y)] <- Inf  # log(Inf)=Inf
  # -------------------------------------------------------------------------- #
  # Adaptive probability
  if((length(b0) != 1) || !is.numeric(b0) || (b0 < 0)) {
    stop("'b0' not valid: it must be greater than or equal to 0.")
  }
  if((length(b1) != 1) || !is.numeric(b1) || (b1 <= 0)) {
    stop("'b1' not valid: it must be greater than 0.")
  }
  prob <- 1 / exp(b0 + b1 * seq(1, nrun))
  uu <- runif(nrun)
  # -------------------------------------------------------------------------- #
  # p constant (c_p)
  if(is.null(p_constant)) {
    p_constant <- 10 * exp(1) * log(p) / p
  } else if((length(p_constant) != 1) || !is.numeric(p_constant) || (p_constant <= 0) || (p_constant >= 1)) {
    stop("'p_constant' not valid: it must be a number in (0,1).")
  }
  # -------------------------------------------------------------------------- #
  # check fixed parameters
  if((length(sd_gammaT) != 1) || !is.numeric(sd_gammaT) || (sd_gammaT <= 0)) {
    stop("'sd_gammaT' not valid: it must be greater than 0.")
  }
  if((length(sd_beta) != 1) || !is.numeric(sd_beta) || (sd_beta <= 0)) {
    stop("'sd_beta' not valid: it must be greater than 0.")
  }
  if((length(sd_gammaB) != 1) || !is.numeric(sd_gammaB) || (sd_gammaB <= 0)) {
    stop("'sd_gammaB' not valid: it must be greater than 0.")
  }
  if((length(a_theta) != 1) || !is.numeric(a_theta) || (a_theta <= 0)) {
    stop("'a_theta' not valid: it must be greater than 0.")
  }
  if((length(b_theta) != 1) || !is.numeric(b_theta) || (b_theta <= 0)) {
    stop("'b_theta' not valid: it must be greater than 0.")
  }
  if((length(a_sigma) != 1) || !is.numeric(a_sigma) || (a_sigma <= 0)) {
    stop("'a_sigma' not valid: it must be greater than 0.")
  }
  if((length(b_sigma) != 1) || !is.numeric(b_sigma) || (b_sigma <= 0)) {
    stop("'b_sigma' not valid: it must be greater than 0.")
  }
  if((length(alpha) != 1) || !is.numeric(alpha) || (alpha < 0)) {
    stop("'alpha' not valid: it must be greater than or equal to 0.")
  }
  # -------------------------------------------------------------------------- #
  # Initialize sigma^-2
  ps <- rgamma(p, a_sigma, b_sigma)
  # Initialize parameters related to the covariates, if x exists
  if(!xnull) {
    Beta <- matrix(rnorm(d * p, 0, sd_beta), nrow = p, ncol = d)    # mean coeff of the data
    GammaT <- matrix(rnorm(qT * d), nrow = qT, ncol = d)            # x effects on Beta coeff
    # precision of GammaT and Beta
    prec_gammaT  <- 1 / (sd_gammaT)  ^ 2
    prec_beta <- 1 / (sd_beta) ^ 2
  } else {
    # not used but necessary to call the C++ function
    Beta <- GammaT <- matrix(1, nrow = 1, ncol = 1)
    prec_gammaT <- prec_beta <- 0.1
  }
  # Initialize lambda star (pxq)
  Lambda_star <- matrix(rnorm(p * k), nrow = p, ncol = k)  # loading matrix
  # Initialize eta (nxk)
  eta <- matrix(rnorm(n * k), nrow = n, ncol = k)          # latent factors
  # Initialize GammaB (qxk) and pred (pxk)
  GammaB <- matrix(rnorm(qB * k), nrow = qB, ncol = k)     # traits effect on local shrinkage
  pred <- wB %*% GammaB                                    # local shrinkage coefficients
  logit <- plogis(pred)
  # Initialize Phi pxk
  Phi <- matrix(rbinom(p * k, size = 1, prob = p_constant), nrow = p, ncol = k)
  # Initialize pi_h, h = 1, ..., k
  v <- c(rbeta(k - 1, shape1 = 1, shape2 = alpha), 1)
  w <- v * c(1, cumprod(1 - v[-k]))  # product up to  l - 1
  d <- rep(k - 1, k)                 # augmented data
  rho <- rep(1, k)                   # preallocation for Bernoulli
  # Initialize the precision matrix of lambda star
  Plam <- diag(rgamma(k, a_theta, b_theta))
  # Compute Lambda (pxk)
  Lambda <- t(t(Lambda_star) * sqrt(rho)) * sqrt(Phi)
  # -------------------------------------------------------------------------- #
  # Allocate output object memory
  valid_outputs <- c("beta",        # coefSamples          : pxd
                     "gammaT",      # GammaTCoefSamples    : qTxd
                     "gammaB",      # shrinkCoefSamples    : qBxk
                     "eta",         # etaval               : nxk
                     "lambda",      # loadSamples          : pxk
                     "sigmacol"     # sigmacol (1/sigma^2) : p
  )
  if("all" %in% output) {
    output <- valid_outputs
  } else {
    output <- intersect(output, valid_outputs)
  }
  # -------------------------------------------------------------------------- #
  out <- list("numFactors" = NA)
  if(!xnull) {
    if("beta" %in% output) out["beta"] <- NA
    if("gammaT" %in% output) out["gammaT"] <- NA
  }
  if("gammaB" %in% output) out["gammaB"] <- NA
  if("eta" %in% output) out["eta"] <- NA
  if("lambda" %in% output) out["lambda"] <- NA
  if("sigmacol" %in% output) out["sigmacol"] <- NA
  # -------------------------------------------------------------------------- #
  # start time
  t0 <- proc.time()
  # -------------------------------------------------------------------------- #
  # ADAPTIVE GIBBS SAMPLING
  # -------------------------------------------------------------------------- #
  out <- Rcpp_cosin(alpha, a_sigma, a_y, a_yp1, a_theta,
                    GammaB, b_sigma, b0, b1, burn, GammaT, b_theta,
                    d,
                    eta,
                    kmax, kstar,
                    Lambda, Lambda_star, logit,
                    Beta,
                    nrun,
                    out,
                    Phi, Plam, prec_gammaT, prec_beta, pred, prob, ps, p_constant,
                    rho,
                    sd_gammaB, sp, start_adapt,
                    thin,
                    uu,
                    v, verbose,
                    w, x, xnull,
                    wT, wB)
  # -------------------------------------------------------------------------- #
  for (it in 1:sp) {
    if("sigmacol" %in% output) out$sigmacol[[it]] <- c(out$sigmacol[[it]])
  }
  out[["numFactors"]] <- c(out[["numFactors"]])
  out[["time"]] <- (proc.time() - t0)[1]
  out[["y"]] <- y                 # data                       : nxp
  if(!xnull) out[["x"]] <- x      # covariates                 : nxd
  out[["wT"]] <- wT               # technical meta-covariates  : pxqT
  out[["wB"]]  <- wB              # biological meta-covariates : pxqB
  out[["hyperparameters"]] <- list(alpha = alpha, a_theta = a_theta,
                                   b_theta = b_theta, sd_gammaT = sd_gammaT,
                                   sd_beta = sd_beta, sd_gammaB = sd_gammaB,
                                   a_sigma = a_sigma, b_sigma = b_sigma, p_constant = p_constant,
                                   xFormula = xFormula, wTformula = wTformula,
                                   wBformula = wBformula, y_max = y_max)
  # -------------------------------------------------------------------------- #
  return(out)
  # -------------------------------------------------------------------------- #
}
