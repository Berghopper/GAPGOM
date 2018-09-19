# Main function is Prediction Function;
#
#                       Input: 1- GeneID - Ensembl ID of the gene that you want to predict, e.g. ENSG00000228630
#                              2- Onto - Ontology type, currently two options; MF(Molecular Function) or BP(Biological Process)
#                              3- Method - Currently five options; Pearson, Spearman, Fisher, Sobolev, combine
#
#                       Output: the output is a list of ontology terms, ordered with respect to FDR values
#                              1- GOID - Gene Ontology ID
#                              2- Ontology - Ontology type (MF or BP)
#                              3- FDR - False Positive Rate
#                              4- Term - description of GOID
#



#' UGenePred - prediction_function()
#'
#' Predicts functions of un-annotated genes based on Gene Ontology and correlated expression patterns.
#'
#' @param gene_id gene to be compared to the other GO terms.
#' @param ontology desired ontology to use for prediction
#' @param expression_data lncrna2function expression data matrix
#' @param ensembl_to_go_id_conversion_df matrix with corresponding GOIDs and EnsemblIDs.
#' @param start_expression_col should be the starting column containing expression values.
#' @param end_expression_col should be the last column containing expression values.
#' @param enrichment_cutoff default is 250
#' @param method default is "combine"
#' @return The resulting matrix with prediction of similar GO terms.
#' @examples UGenePred::prediction_function("ENSG00000228630", "BP", ExpressionData, EnsemblID2GOID, 4, ncol(ExpressionData))
#'
#' @export
prediction_function <- function(gene_id, ontology, expression_data=UGenePred::expression_data, ensembl_to_go_id_conversion_df=UGenePred::ensembl_id_to_go_id, start_expression_col=4, end_expression_col=22, enrichment_cutoff = 250, method = "combine") {
  # prepare the data with some special operations/vars that are needed later
  expression_data_pc <- expression_data[(expression_data$GeneType == "protein_coding"), ]
  ensembl_id_pc <- expression_data_pc$GeneID

  # Target expression data where gene id matches
  target_expression_data <- expression_data[(expression_data[, 1] == gene_id), ]
  target_expression_data <- base::as.numeric(target_expression_data[1, c(4:end_expression_col)])

  # these functions calculate score between target expression of target gene vs the rest of the desired protein coding genes.  this is either a correlation metric or a geometrical
  # metric.  after score is calculated, these scores are enriched.
  enrichment_result = NULL
  if (method == "Pearson") {
    score <- base::apply(expression_data_pc[, c(start_expression_col:end_expression_col)], 1, function(x) base::abs(stats::cor(base::as.numeric(x), target_expression_data)))
    score_df <- UGenePred::prepare_score_df(ensembl_id_pc, score, gene_id)
    ordered_score_df <- score_df[base::rev(base::order(score_df[, 2])), ]
    enrichment_result <- UGenePred::enrichment_analysis(ordered_score_df, ontology, ensembl_id_pc=ensembl_id_pc, ensembl_to_go_id_conversion_df=ensembl_to_go_id_conversion_df, enrichment_cutoff=enrichment_cutoff)

  } else if (method == "Spearman") {
    score <- base::apply(expression_data_pc[, c(start_expression_col:end_expression_col)], 1, function(x) base::abs(stats::cor(base::as.numeric(x), target_expression_data, method = "spearman")))
    score_df <- UGenePred::prepare_score_df(ensembl_id_pc, score, gene_id)
    ordered_score_df <- score_df[base::rev(base::order(score_df[, 2])), ]
    enrichment_result <- UGenePred::enrichment_analysis(ordered_score_df, ontology, ensembl_id_pc=ensembl_id_pc, ensembl_to_go_id_conversion_df=ensembl_to_go_id_conversion_df, enrichment_cutoff=enrichment_cutoff)

  } else if (method == "FisherMetric") {
    score <- base::apply(expression_data_pc[, c(start_expression_col:end_expression_col)], 1, function(x) UGenePred::fisher_metric(as.numeric(x), target_expression_data))
    score_df <- UGenePred::prepare_score_df(ensembl_id_pc, score, gene_id)
    ordered_score_df <- score_df[base::order(score_df[, 2]), ]
    enrichment_result <- UGenePred::enrichment_analysis(ordered_score_df, ontology, ensembl_id_pc=ensembl_id_pc, ensembl_to_go_id_conversion_df=ensembl_to_go_id_conversion_df, enrichment_cutoff=enrichment_cutoff)

  } else if (method == "SobolevMetric") {
    score <- base::apply(expression_data_pc[, c(start_expression_col:end_expression_col)], 1, function(x) UGenePred::sobolev_metric(as.numeric(x), target_expression_data))
    score_df <- UGenePred::prepare_score_df(ensembl_id_pc, score, gene_id)
    ordered_score_df <- score_df[base::order(score_df[, 2]), ]
    enrichment_result <- UGenePred::enrichment_analysis(ordered_score_df, ontology, ensembl_id_pc=ensembl_id_pc, ensembl_to_go_id_conversion_df=ensembl_to_go_id_conversion_df, enrichment_cutoff=enrichment_cutoff)

  } else if (method == "combine") {
    score_pearson <- base::apply(expression_data_pc[, c(start_expression_col:end_expression_col)], 1, function(x) base::abs(stats::cor(base::as.numeric(x), target_expression_data)))
    score_spearman <- base::apply(expression_data_pc[, c(start_expression_col:end_expression_col)], 1, function(x) base::abs(stats::cor(base::as.numeric(x), target_expression_data,
                                                                                                                                        method = "spearman")))
    score_fisher <- base::apply(expression_data_pc[, c(start_expression_col:end_expression_col)], 1, function(x) UGenePred::fisher_metric(base::as.numeric(x), target_expression_data))
    score_sobolev <- base::apply(expression_data_pc[, c(start_expression_col:end_expression_col)], 1, function(x) UGenePred::sobolev_metric(base::as.numeric(x), target_expression_data))

    score_df <- base::data.frame(ensembl_id_pc, score_pearson, score_spearman, score_sobolev, score_fisher)
    score_df <- stats::na.omit(score_df)
    # filter gene ID's (Select everything except the chosen gene).
    score_df <- score_df[(score_df[, 1] != gene_id), ]

    # order each score according to their respective method.
    ordered_score_df_pearson <- score_df[base::rev(base::order(score_df$score_pearson)), ]
    ordered_score_df_spearman <- score_df[base::rev(base::order(score_df$score_spearman)), ]
    ordered_score_df_sobolev <- score_df[base::order(score_df$score_sobolev), ]
    ordered_score_df_fisher <- score_df[base::order(score_df$score_fisher), ]

    # Run enrichment for each method
    enrichment_pearson <- UGenePred::enrichment_analysis(ordered_score_df_pearson, ontology, ensembl_id_pc=ensembl_id_pc, ensembl_to_go_id_conversion_df=ensembl_to_go_id_conversion_df, enrichment_cutoff=enrichment_cutoff)
    enrichment_spearman <- UGenePred::enrichment_analysis(ordered_score_df_spearman, ontology, ensembl_id_pc=ensembl_id_pc, ensembl_to_go_id_conversion_df=ensembl_to_go_id_conversion_df, enrichment_cutoff=enrichment_cutoff)
    enrichment_sobolev <- UGenePred::enrichment_analysis(ordered_score_df_sobolev, ontology, ensembl_id_pc=ensembl_id_pc, ensembl_to_go_id_conversion_df=ensembl_to_go_id_conversion_df, enrichment_cutoff=enrichment_cutoff)
    enrichment_fisher <- UGenePred::enrichment_analysis(ordered_score_df_fisher, ontology, ensembl_id_pc=ensembl_id_pc, ensembl_to_go_id_conversion_df=ensembl_to_go_id_conversion_df, enrichment_cutoff=enrichment_cutoff)

    # before binding, add column for showing what method result is originated from.
    enrichment_pearson[, "used_method"] <- base::rep("Pearson", base::nrow(enrichment_pearson))
    enrichment_spearman[, "used_method"] <- base::rep("Spearman", base::nrow(enrichment_spearman))
    enrichment_sobolev[, "used_method"] <- base::rep("Sobolev", base::nrow(enrichment_sobolev))
    enrichment_fisher[, "used_method"] <- base::rep("Fisher", base::nrow(enrichment_fisher))
    # bind all results by row
    combined_enrichment <- base::rbind(enrichment_pearson, enrichment_spearman, enrichment_sobolev, enrichment_fisher)

    # and keep the rows with the lowest corrected P-values.
    combined_enrichment <- plyr::ddply(combined_enrichment, plyr::.(GOID), function(x) x[base::which.min(x$FDR), ])
    # Order p-values on lowest > bigest
    enrichment_result <- combined_enrichment[base::order(base::as.numeric(combined_enrichment$FDR)), ]
  } else {
    cat("Selected method; \"", method, "\" does not exist! Refer to the help page for options.")
    return(NULL)
  }
  ###

  if (nrow(enrichment_result) > 0) {
    # number the rownames and return the enrichment results.
    base::rownames(enrichment_result) <- c(1:base::nrow(enrichment_result))
    return(enrichment_result)
  } else {
    print("Could not find any similar genes!")
  }
}
