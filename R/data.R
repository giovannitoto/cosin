#' @title Cell-specific covariates
#' @description A \eqn{n\times d} matrix of cell-specific covariates.
#' @format A data frame with 20 rows and 4 columns:
#' \itemize{
#'   \item \code{unaligned}
#'   \item \code{aligned_unmapped}
#'   \item \code{number_of_genes}
#'   \item \code{total_count_per_cell}
#' }
# @details
"x"

#' @title Gene-specific technical meta-covariates
#' @description A \eqn{p\times q_T} matrix of gene-specific technical meta-covariates.
#' @format A data frame with 100 rows and 2 columns:
#' \itemize{
#'   \item \code{length}
#'   \item \code{gc_content}
#' }
# @details
"wT"

#' @title Gene-specific biological meta-covariates
#' @description A \eqn{p\times q_B} matrix \eqn{x} of gene-specific biological meta-covariates.
#' @format A data frame with 100 rows and 3 columns:
#' \itemize{
#'   \item \code{Alzheimer_disease}
#'   \item \code{Parkinson_disease}
#'   \item \code{Huntington_disease}
#' }
# @details
"wB"

#' @title Matrix of gene expression
#' @description A \eqn{n\times p} matrix \eqn{y} of gene expression containing the counts
#' for \eqn{p} genes measured in on \eqn{n} cells.
#' @format A matrix with 20 rows and 100 columns.
# @details
"y"
