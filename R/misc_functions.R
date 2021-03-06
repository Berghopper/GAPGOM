#' GAPGOM internal - Miscellaneous functions
#'
#' These functions are misc. functions that are mostly used in applies. They
#' are defined here probably because they are used often.
#'
#' @param x vector of expression values
#' @param args internal arguments, used for target expression datas
#'
#' @section Notes:
#' These functions are internal functions and should not be called by the user.
#'
#' @return output is different on a case-to-case basis
#' 
#' @name misc_functions
#' @keywords internal
NULL

#' @rdname misc_functions
misc_pearson <- function(x, args) {
  return(abs(cor(as.numeric(x), as.numeric(args$target_expression_data))))}

#' @rdname misc_functions
misc_sobolev <- function(x, args) {
  return(sobolev_metric(as.numeric(x), as.numeric(args$target_expression_data)))}

#' @rdname misc_functions
misc_fisher <- function(x, args) {
  return(fisher_metric(as.numeric(x), as.numeric(args$target_expression_data)))}

#' @rdname misc_functions
misc_kendall <- function(x, args) {
  return(abs(cor(as.numeric(x), as.numeric(args$target_expression_data),
                 method = "kendall")))}

#' @rdname misc_functions
misc_spearman <- function(x, args) {
  return(abs(cor(as.numeric(x), as.numeric(args$target_expression_data),
                 method = "spearman")))}
