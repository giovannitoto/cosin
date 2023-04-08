# ---------------------------------------------------------------------------- #

#' Contributions
#'
#' Compute the sorted contributions \eqn{C^*_1,\ldots,C^*_{k_{max}}} for each iteration.
#'
#' @param out_MCMC A list containing the results of an Adaptive Gibbs Sampler obtained using the function \code{\link{cosin}}.
#' @param reference An integer specifying the reference iteration used to sort the contributions of each iteration; if \code{reference=NULL}, the contributions are sorted using the Frobenius norm. Default is \code{NULL}.
#' @param chains Logical: if \code{TRUE}, the function returns also the draws from the posterior distribution of \eqn{C^*_1,\ldots,C^*_{k_{max}}}.
#' @param verbose Logical: if \code{TRUE}, print the number of active factors every 50 iterations. Default is \code{TRUE}.
#'
#' @return A list with the following elements:
#' \itemize{
#' \item \code{C1}: a list of matrices containing the posterior means of \eqn{C^*_1,\ldots,C^*_{k_{max}}}.
#' \item \code{C2}: a list of matrices containing the posterior means of the square of the elements of \eqn{C^*_1,\ldots,C^*_{k_{max}}}.
#' \item \code{reference}: index of the iteration used as reference.
#' \item \code{chains}: a list containing draws from the posterior distribution of \eqn{C^*_1,\ldots,C^*_{k_{max}}}.
#' }
#'
#'
#' @seealso This function is applied to an output of \code{\link{cosin}}.
#'
#' @export
contributions <- function(out_MCMC, reference = NULL, chains = FALSE, verbose = TRUE) {
  # -------------------------------------------------------------------------- #
  # check whether all necessary variables are available in out_MCMC
  required_variables <- c("numFactors", "eta", "lambda")
  if(check_list(out_MCMC, required_variables) == FALSE) {
    stop(paste(paste(c(required_variables), collapse = ", "), "must be stored in the MCMC output."))
  }
  # -------------------------------------------------------------------------- #
  sp <- length(out_MCMC$numFactors)
  n <- nrow(out_MCMC$eta[[1]])
  p <- nrow(out_MCMC$lambda[[1]])
  k_max <- max(sapply(out_MCMC$lambda, ncol))
  # -------------------------------------------------------------------------- #
  if(!is.null(reference)) {
    if((length(reference) != 1) || !is.numeric(reference) || (reference != round(reference)) || reference > sp) {
      stop("'reference' not valid: it must be NULL or an integer between 1 and the number of iterations.")
    }
    # number of factors of reference iteration
    k_ref <- ncol(out_MCMC$lambda[[reference]])
    # compute contributions NOT in order of reference iteration
    C_list <- list()
    for (h in 1:k_ref) {
      C_list[[paste0("C", h)]] <- tcrossprod(out_MCMC$eta[[reference]][, h], out_MCMC$lambda[[reference]][, h])
    }
    # compute order of the contributions using Frobenius norm
    C_order <- order(sapply(C_list, function(C) norm(C, type = "F")), decreasing = TRUE)
    # order contributions of reference iteration
    C_ref <- list()
    for (h in 1:k_ref) {
      C_ref[[paste0("C", h)]] <- C_list[[paste0("C", C_order[h])]]
    }
  }
  # -------------------------------------------------------------------------- #
  if((length(verbose) != 1) || !is.logical(verbose)) {
    stop("'verbose' not valid: it must be 'TRUE' or 'FALSE'.")
  }
  # -------------------------------------------------------------------------- #
  if((length(chains) != 1) || !is.logical(chains)) {
    stop("'chains' not valid: it must be 'TRUE' or 'FALSE'.")
  }
  if(chains) {
    chains <- list()
    for (h in 1:k_max) {
      chains[[paste0("C", h)]] <- vector("list", length = sp)
    }
    save <- TRUE
  } else {
    save <- FALSE
  }
  # -------------------------------------------------------------------------- #
  C1 <- C2 <- list()
  for (h in 1:k_max) {
    C1[[paste0("C", h)]] <- C2[[paste0("C", h)]] <- matrix(0, nrow = n, ncol = p)
  }
  # -------------------------------------------------------------------------- #
  for (it in 1:sp) {
    # number of factors
    k <- ncol(out_MCMC$lambda[[it]])
    # compute contributions NOT in order
    C_list <- list()
    for (h in 1:k) {
      C_list[[paste0("C", h)]] <- tcrossprod(out_MCMC$eta[[it]][, h], out_MCMC$lambda[[it]][, h])
    }
    # add matrices of zeros if k is less than the maximum k observed
    if (k < k_max) {
      for (h in (k+1):k_max) C_list[[paste0("C", h)]] <- matrix(0, nrow = n, ncol = p)
    }
    # compute order using Frobenius norm or reference iteration
    if(is.null(reference)) {
      # compute order of the contributions using Frobenius norm
      C_order <- order(sapply(C_list, function(C) norm(C, type = "F")), decreasing = TRUE)
      for (h in 1:k_max) {
        # update C1 and C2
        C1[[paste0("C", h)]] <- C1[[paste0("C", h)]] + C_list[[paste0("C", C_order[h])]] / sp
        C2[[paste0("C", h)]] <- C2[[paste0("C", h)]] + C_list[[paste0("C", C_order[h])]]^2 / sp
        # save contributions
        if(save) chains[[paste0("C", h)]][[it]] <- C_list[[paste0("C", C_order[h])]]
      }
    } else {
      # compute order of the contributions using reference iteration
      for (h in 1:k_ref) {
        min_idx <- names(which.min(sapply(C_list, function(C) norm(C - C_ref[[paste0("C", h)]], type = "F"))))
        # update mean and second moment
        C1[[paste0("C", h)]] <- C1[[paste0("C", h)]] + C_list[[min_idx]] / sp
        C2[[paste0("C", h)]] <- C2[[paste0("C", h)]] + C_list[[min_idx]]^2 / sp
        # save contribution
        if(save) chains[[paste0("C", h)]][[it]] <- C_list[[min_idx]]
        C_list[[min_idx]] <- NULL
      }
      if(k_ref < k_max) {
        # the remaining contributions are sorted using Frobenius norm
        for (h in (k_ref+1):k_max) {
          min_idx <- names(which.min(sapply(C_list, function(C) norm(C, type = "F"))))
          # update C1 and C2
          C1[[paste0("C", h)]] <- C1[[paste0("C", h)]] + C_list[[min_idx]] / sp
          C2[[paste0("C", h)]] <- C2[[paste0("C", h)]] + C_list[[min_idx]]^2 / sp
          # save contribution
          if(save) chains[[paste0("C", h)]][[it]] <- C_list[[paste0("C", min_idx)]]
          C_list[[paste0("C", min_idx)]] <- NULL
        }
      }
    }
    if(verbose == TRUE && it %% 100 == 0) cat(it, ":", out_MCMC$numFactors[it], "active factors\n")
  }
  # -------------------------------------------------------------------------- #
  output <- list("C1" = C1, "C2" = C2, "reference" = reference)
  if(save) output[["chains"]] <- chains
  # -------------------------------------------------------------------------- #
  return(output)
  # -------------------------------------------------------------------------- #
}

# ---------------------------------------------------------------------------- #
