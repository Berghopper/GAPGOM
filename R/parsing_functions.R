#' UGenePred internal - prepare_score_df()
#'
#' This function is an internal function and should not be called by the user.
#'
#' @section Notes:
#' internal function used in prediction_function()
#'
#' @param pc_ensmlb_id Ensmlb ID's for protein coding expression data.
#' @param score score vector
#' @param gene_id Gene ID originally selected in the prediction function.
#' @return
prepare_score_df <- function(pc_ensmbl_id, score, gene_id) {
    # score_df is the 'score' (correlation/metric) dataframe against target gene.  combine the dataframe if combine is selected, otherwise just add regular score.
    score_df <- base::data.frame(pc_ensmbl_id, score)
    score_df <- stats::na.omit(score_df)
    # filter gene ID's (Select everything except the chosen gene).
    score_df <- score_df[(score_df[, 1] != gene_id), ]
    return(score_df)
}

# convert an ext ID to a term ID dat = EG_ (Extracted genes with correct gene ontology and narrowed selection of top genes.)  List_Top_Genes (List of top 250 genes)
ext_id_to_term_id <- function(dat, List_Top_Genes) {
    # add 1 for each gene that is in list top genes.
    return(plyr::ddply(dat, plyr::.(go_id), function(x) base::nrow(x[(x$ensembl_gene_id %in% List_Top_Genes), ])))
}
