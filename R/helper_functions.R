# ---------------------------------------------------------------------------- #

#' Check if all required variables are a list
#'
#' @param list_to_check A list.
#' @param required_variables A vector containing the names of the variables which must be in \code{list_to_check}.
#'
#' @return It returns \code{TRUE} if \code{list_to_check} contains all the required variables specified in \code{required_variables}.
check_list <- function(list_to_check, required_variables) {
  if(all(required_variables %in% names(list_to_check))) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

# ---------------------------------------------------------------------------- #

#' Inversion of omega
#'
#' This function is called in 'posterior_mean' and 'lposterior_function'.
#'
#' @param lambda A pxk matrix of factorial weights.
#' @param sigmai A p-dimensional vector of \code{1/sigma^2}.
#'
#' @return A pxp matrix.
omega_inversion <- function(lambda, sigmai) {
  p <- dim(lambda)[1]
  k <- dim(lambda)[2]
  # Woodbury matrix identity
  lambda_sigmai <- crossprod(lambda, diag(sigmai))
  capac <- diag(k) + lambda_sigmai %*% lambda
  omega_inv <- diag(sigmai) - diag(sigmai) %*% lambda %*% solve(capac) %*% lambda_sigmai
  return(omega_inv)
}

# ---------------------------------------------------------------------------- #

#' Rescale one or more parameters
#'
#' Add columns to one or more parameters whose number of columns depends on the number of factors so that the parameters of all MCMC iterations have the same number of columns. This number can be equal to the maximum number of active factors or the maximum number of factors, active and inactive, observed during the Adaptive Gibbs Sampler.
#'
#' @param out_MCMC A list containing the results of an Adaptive Gibbs Sampler obtained using the function \code{\link{cosin}}.
#' @param parameters A vector of strings specifying the parameters for which you want to compute the posterior mean. The possible valid strings are \code{gammaB}, \code{eta} and \code{lambda}. Default is \code{all}, which is equivalent to writing \code{c("gammaB", "eta", "lambda")}.
#' @param columns A string specifying whether to consider only active factors (\code{kstar}) or both active and inactive factors (\code{k}). Default is \code{k}.
#' @param samples A vector of integers specifying the iterations to consider. By default, all iterations are considered.
#'
#' @return A list containing the MCMC chains of the parameters specified in \code{parameters}.
#'
#' @seealso This function is applied to an output of \code{\link{cosin}}.
rescale_parameter <- function(out_MCMC, parameters = "all", columns = "k", samples = NULL) {
  # check the value of 'columns'
  if(!(columns %in% c("k", "kstar"))) {
    stop("'columns' not valid: it must be 'k' or 'kstar'.")
  }
  # remove invalid strings from 'parameters'
  valid_parameters <- c("gammaB", "eta", "lambda")
  if("all" %in% parameters) {
    parameters <- valid_parameters
  } else {
    parameters <- intersect(parameters, valid_parameters)
  }
  # stop if no valid string is provided as input
  if(length(parameters) == 0) {
    stop("No valid parameter provided as input.")
  }
  # check whether all necessary variables are available in out_MCMC
  required_variables <- c("numFactors", parameters)
  if(check_list(out_MCMC, required_variables) == FALSE) {
    stop(paste(paste(required_variables, collapse = ", "), "must be stored in the MCMC output."))
  }
  # number of iterations
  t <- length(out_MCMC$numFactors)
  # check if the indices of selected iterations are correct
  if (is.null(samples)) {
    samples <- 1:t
  } else {
    if(min(samples) < 1 | max(samples) > t) {
      stop(paste("'samples' not valid: it must be a vector of integers between 1 and ", t, ".", sep=""))
    }
    if(any(sapply(samples, function(x) x != round(x)))) {
      stop("'samples' not valid: it must be a vector of integers.")
    }
  }
  # -------------------------------------------------------------------------- #
  if("gammaB" %in% parameters) q <- nrow(out_MCMC$gammaB[[1]])
  if("eta" %in% parameters) n <- nrow(out_MCMC$eta[[1]])
  if("lambda" %in% parameters) p <- nrow(out_MCMC$lambda[[1]])
  # prepare output
  output <- list()
  for (par in parameters) output[[par]] <- list()
  # max number of factors/columns
  if(columns == "kstar") {
    K <- max(out_MCMC$numFactors)
  } else {
    K <- max(sapply(out_MCMC[[names(output)[1]]], ncol))
  }
  # complete the parameters with the correct number of factors for each iteration
  for (i in samples) {
    # ------------------------------------------------------------------------ #
    if(columns == "kstar") {
      if("gammaB" %in% parameters) out_MCMC$gammaB[[i]] <- out_MCMC$gammaB[[i]][, 1:out_MCMC$numFactors[i]]
      if("eta" %in% parameters) out_MCMC$eta[[i]] <- out_MCMC$eta[[i]][, 1:out_MCMC$numFactors[i]]
      if("lambda" %in% parameters) out_MCMC$lambda[[i]] <- out_MCMC$lambda[[i]][, 1:out_MCMC$numFactors[i]]
    }
    # ------------------------------------------------------------------------ #
    kl <- ncol(out_MCMC[[names(output)[1]]][[i]])
    if(kl > K) {
      # if the number of factors is greater than K, keep only the first 'K' factors/columns
      if("gammaB" %in% parameters) output$gammaB[[i]] <- out_MCMC$gammaB[[i]][, 1:K]
      if("eta" %in% parameters) output$eta[[i]] <- out_MCMC$eta[[i]][, 1:K]
      if("lambda" %in% parameters) output$lambda[[i]] <- out_MCMC$lambda[[i]][, 1:K]
    } else {
      # if the number of factors is less than K, add columns of 0s until you have 'K' factors/columns
      if("gammaB" %in% parameters) {
        vec <- matrix(rep(0, q * (K - kl)), nrow = q, ncol = K - kl)
        output$gammaB[[i]] <- cbind(out_MCMC$gammaB[[i]], vec)         # qxkmax
      }
      if("eta" %in% parameters) {
        vec <- matrix(rep(0, n * (K - kl)), nrow = n, ncol = K - kl)
        output$eta[[i]] <- cbind(out_MCMC$eta[[i]], vec)           # nxkmax
      }
      if("lambda" %in% parameters) {
        vec <- matrix(rep(0, p * (K - kl)), nrow = p, ncol = K - kl)
        output$lambda[[i]] <- cbind(out_MCMC$lambda[[i]], vec)     # pxkmax
      }
    }
    # ------------------------------------------------------------------------ #
  }
  return(output)
}

# ---------------------------------------------------------------------------- #
