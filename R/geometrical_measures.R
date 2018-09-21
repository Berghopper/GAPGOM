#' UGenePred - Geometrical Similarity Measures
#'
#' These geometric functions will calculate the distance between two random
#' sets of variables. In our usecase this would be an expression values (FPKM).
#' Given two gene expression value lists, semantic similarity will be
#' calculated between them. This is then returned.
#'
#' Fisher metric is based on the G. Lebanon et al. implementation. [1]
#'
#' The sobolev metrix is based on the T. Villmann et al. implementation. [2]
#'
#' @section Notes:
#' Both sobolev_metric() and fisher_metric() should have exactly the same
#' input and output types.
#'
#' @param x gene expression value vector (numeric) of reference/novel gene
#' that needs to be compared to the dataset.
#' @param y also a gene expression value vector, the idea here is that indexes
#' between x and y correspond for same or similar tissues for accurate
#' similarity. x and y may also be swapped.
#' @return score of similarity between the two vectors (type=double)
#' @references [1]. Villmann T: \strong{Sobolev metrixs for learning of functional
#' data - mathematical and theoretical aspects.} In:
#' \emph{Machine Learning Reports}.
#' Edited by Villmann T, Schleif F-m, vol. 1. Leipzig, Germany: Medical
#' Department, University of Leipzig; 2007: 1-13.
#' @references [2]. Lebanon G: \strong{Learning riemannian metrics.} In:
#' \emph{Proceedings of the Nineteenth conference on Uncertainty in Artificial
#' Intelligence; Acapulco Mexico.} Morgan Kaufmann Publishers Inc. 2003: 362-369
#' @name geometrical_measures
NULL


#' @rdname geometrical_measures
#' @export
sobolev_metric <- function(x, y) {
    x1 <- x^2/(base::sum(x^2))
    y1 <- y^2/(base::sum(y^2))
    z1 <- x1 - y1
    FT <- stats::fft(z1)
    w <- 2 * pi * (1:base::length(FT))/(base::length(FT))
    s <- base::sum((1 + w) * base::abs(FT)^2)^(1/2)
    return(s)
}

#' @rdname geometrical_measures
#' @export
fisher_metric <- function(x, y) {
    x1 <- x^2/(base::sum(x^2))
    y1 <- y^2/(base::sum(y^2))
    t <- x1 * y1
    s <- base::acos(base::sum(base::sqrt(t)))
    return(s)
}
