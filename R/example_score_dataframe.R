#' GAPGOM - example similarity score data
#'
#' A dataframe containing similarity scores towards another gene of choice.
#' The gene that was compared to the rest in this case is "ENSG00000228630".
#' Similarity scores for this dataframe were obtained using the Fisher method.
#'
#' @format A data frame with 20183 rows and 2 variables/columns:
#' \describe{
#'   \item{pc_ensembl_id}{EnsemblID of protein coding gene.}
#'   \item{GeneType}{calculated similarity score with measure of choice, for
#'   this data, Fisher was used.}
#' }
#'
#' @source Intermediate score_df output of prediction_function() in internal
#' predict_fisher() method
"example_score_dataframe"
