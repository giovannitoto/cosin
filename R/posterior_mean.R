# ---------------------------------------------------------------------------- #

#' Posterior mean
#'
#' Compute the posterior mean of one or more parameters.
#'
#' @param out_MCMC A list containing the results of an Adaptive Gibbs Sampler obtained using the function \code{\link{cosin}}.
#' @param parameters A vector containing the names of the parameters for which you want to compute the posterior mean. The possible valid strings are \code{"beta"}, \code{"gammaT"}, \code{"sigmacol"}, \code{"omega"}, \code{"omega_inv"} and \code{"omega_pcorr"}. Default is \code{"all"}, which is equivalent to writing \code{c("beta", "gammaT", "sigmacol", "omega", "omega_inv", "omega_pcorr")}.
#' @param columns A string specifying whether to consider only active factors (\code{"kstar"}) or both active and inactive factors (\code{"k"}) for the computation of \code{omega}, \code{omega_inv} and \code{"omega_pcorr"}. Default is \code{"k"}.
#'
#' @return A list containing the posterior means of the parameters specified in \code{parameters}.
#'
#' @seealso This function is applied to an output of \code{\link{cosin}}.
#'
#' @export
posterior_mean <- function(out_MCMC, parameters = "all", columns = "k") {
  # remove invalid strings from 'parameters'
  valid_parameters <- c("beta", "gammaT", "sigmacol", "omega", "omega_inv", "omega_pcorr")
  if("all" %in% parameters) {
    parameters <- valid_parameters
  } else {
    parameters <- intersect(parameters, valid_parameters)
  }
  # stop if no valid string is provided as input
  if(length(parameters) == 0) {
    stop("No valid parameter provided as input.")
  }
  # prepare output
  output <- list()
  # check whether all necessary variables are available in out_MCMC
  for (par in intersect(parameters, valid_parameters[1:3])) {
    if(par %in% names(out_MCMC)) {
      # compute the posterior mean of the parameter 'par'
      output[[par]] <- Reduce("+", out_MCMC[[par]]) / length(out_MCMC[[par]])
    } else {
      stop(paste("'", par, "' not stored in the MCMC output.", sep=""))
    }
  }
  # compute the posterior mean of 'omega' and/or 'omega_inv'
  if (any(c("omega", "omega_inv", "omega_pcorr") %in% parameters)) {
    # check the value of 'columns'
    if((length(columns) != 1) || !(columns %in% c("k", "kstar"))) {
      stop("'columns' not valid: it must be 'k' or 'kstar'.")
    }
    if(check_list(out_MCMC, c("lambda", "numFactors", "sigmacol")) == FALSE) {
      stop("'lambda', 'numFactors' and 'sigmacol' must be stored in the MCMC output.")
    } else {
      # number of iterations
      t <- length(out_MCMC$numFactors)
      p <- nrow(out_MCMC$lambda[[1]])
      out_MCMC$lambda <- rescale_parameter(out_MCMC = out_MCMC, parameters = "lambda", columns = columns)$lambda
      if ("omega" %in% parameters) {
        output[["omega"]] <- matrix(0, nrow = p, ncol = p)
        for (i in 1:t) {
          output[["omega"]] <- output[["omega"]] + tcrossprod(out_MCMC$lambda[[i]]) + diag(1 / out_MCMC$sigmacol[[i]])
        }
        output[["omega"]] <- output[["omega"]] / t
      }
      if ("omega_inv" %in% parameters) {
        output[["omega_inv"]] <- matrix(0, nrow = p, ncol = p)
        for (i in 1:t) {
          output[["omega_inv"]] <- output[["omega_inv"]] + omega_inversion(out_MCMC$lambda[[i]], out_MCMC$sigmacol[[i]])
        }
        output[["omega_inv"]] <- output[["omega_inv"]] / t
      }
      if("omega_pcorr" %in% parameters) {
        output[["omega_pcorr"]] <- matrix(0, nrow = p, ncol = p)
        for (i in 1:t) {
          omega <- tcrossprod(out_MCMC$lambda[[i]]) + diag(1 / out_MCMC$sigmacol[[i]])
          output[["omega_pcorr"]] <- output[["omega_pcorr"]] + stats::cov2cor(solve(omega))
        }
        output[["omega_pcorr"]] <- output[["omega_pcorr"]] / t
      }
    }
  }
  # return a list
  return(output)
}

# ---------------------------------------------------------------------------- #
