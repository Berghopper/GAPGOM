#' UGenePred internal - prepare_score_df()
#'
#' This function is an internal function and should not be called by the user.
#'
#' Prepared score dataframe for the enrichment analysis.
#'
#' @section Notes:
#' internal function used in prediction_function()
#'
#' @param pc_ensmlb_id Ensmlb ID's for protein coding expression data.
#' @param score score vector
#' @param gene_id Gene ID originally selected in the prediction function.
#' @return The score dataframe with ensmbl ID's
prepare_score_df <- function(pc_ensmbl_id, score, gene_id) {
    # score_df is the 'score' (correlation/metric) dataframe against target gene.  combine the dataframe if combine is selected, otherwise just add regular score.
    score_df <- base::data.frame(pc_ensmbl_id, score)
    score_df <- stats::na.omit(score_df)
    # filter gene ID's (Select everything except the chosen gene).
    score_df <- score_df[(score_df[, 1] != gene_id), ]
    return(score_df)
}

#' UGenePred internal - ext_id_to_term_id()
#'
#' This function is an internal function and should not be called by the user.
#'
#' Convert an ext ID to a term ID dat = EG_ (Extracted genes with correct gene ontology and narrowed selection of top genes.)  list_top_genes (List of top 250 genes)
#'
#' @section Notes:
#' internal function used in enrichment_analysis()
#'
#' @param data extracted_genes with correct onto etc. (matrix)
#' @param list_top_genes The list_top_genes matrix
#' @return return the quantified matrix
ext_id_to_term_id <- function(data, list_top_genes) {
    # add 1 for each gene that is in list top genes.
    return(plyr::ddply(data, plyr::.(go_id), function(x) base::nrow(x[(x$ensembl_gene_id %in% list_top_genes), ])))
}
