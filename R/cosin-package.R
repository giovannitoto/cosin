#' @useDynLib cosin, .registration = TRUE
#' @importFrom Rcpp sourceCpp
#' @importFrom Rcpp evalCpp
NULL

.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Welcome to the jungle of Generalized Infinite Factorization models!")
}
