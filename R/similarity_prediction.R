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
# BUG; Error: $ operator is invalid for atomic vectors
#' @export
prediction_function <- function(gene_id, ontology, method="combine", expression_data, start_expression_col, end_expression_col){
  # prepare the data with some special operations/vars that are needed later
  expression_data_pc <- expression_data[(expression_data$GeneType == "protein_coding"), ]
  ensembl_id_pc <- expression_data_pc$GeneID

  # Target expression data where gene id matches
  target_expression_data <- expression_data[( expression_data[,1] == gene_id ),]
  target_expression_data <- base::as.numeric( target_expression_data[1,c(4:end_expression_col)] )

  # these functions calculate score between target expression of target gene vs
  # the rest of the desired protein coding genes.
  # this is either a correlation metric or a geometrical metric.
  ###
  if(method == "Pearson") {
    score <- base::apply(expression_data_pc[, c(start_expression_col:end_expression_col)], 1, function(x) base::abs(stats::cor(base::as.numeric(x), target_expression_data)))
  }
  if(method == "Spearman") {
    score <- base::apply( expression_data_pc[, c(start_expression_col:end_expression_col)], 1, function(x) base::abs(stats::cor(base::as.numeric(x), target_expression_data, method = "spearman")))
  } else if(method == "FisherMetric") {
    score <- base::apply( expression_data_pc[, c(start_expression_col:end_expression_col)], 1, function(x) UGenePred::SobolevMetric(as.numeric(x), target_expression_data))
  } else if(method == "SobolevMetric") {
    score <- base::apply( expression_data_pc[, c(start_expression_col:end_expression_col)], 1, function(x) UGenePred::FisherMetric(as.numeric(x), target_expression_data))
  } else if(method == "combine") {
    score_pearson  <- base::apply( expression_data_pc[, c(start_expression_col:end_expression_col)], 1, function(x) base::abs(stats::cor(base::as.numeric(x), target_expression_data)))
    score_spearman <- base::apply( expression_data_pc[, c(start_expression_col:end_expression_col)], 1, function(x) base::abs(stats::cor(base::as.numeric(x), target_expression_data, method = "spearman")))
    score_fisher   <- base::apply( expression_data_pc[, c(start_expression_col:end_expression_col)], 1, function(x) UGenePred::FisherMetric(base::as.numeric(x), target_expression_data))
    score_sobolev  <- base::apply( expression_data_pc[, c(start_expression_col:end_expression_col)], 1, function(x) UGenePred::SobolevMetric(base::as.numeric(x), target_expression_data))
  } else {
    cat("Selected method; \"", method,
        "\" does not exist! Refer to the help page for options.")
    return(NULL)
  }
  ###

  # score_df is the "score" (correlation/metric) dataframe against target gene.
  # combine the dataframe if combine is selected, otherwise just add regular
  # score.
  if( method == "combine" ) {
    score_df <- base::data.frame(ensembl_id_pc, score_pearson, score_spearman,
                           score_sobolev, score_fisher)
  } else {
    score_df <- base::data.frame(ensembl_id_pc, score)
  }

  score_df <- stats::na.omit(score_df)
  # filter gene ID's (Select everything except the chosen gene).
  score_df <- score_df[(score_df[,1]!=gene_id),]

  # run the enrichments analysis'.
  if(method == "Pearson" | method == "Spearman") {
    enrichment_result <- UGenePred::enrichment_analysis( score_df[base::rev(base::order(score_df[,2])), ], ontology)
  } else if(method == "SobolevMetric" | method == "FisherMetric") {
    enrichment_result <- UGenePred::enrichment_analysis( score_df[base::order(score_df[,2]), ], ontology)
  } else if(method == "combine") {
    # Run enrichment for each method
    enrichment_pearson  <- UGenePred::enrichment_analysis(score_df[base::rev(base::order(score_df$score_pearson)), ], ontology)
    enrichment_spearman <- UGenePred::enrichment_analysis(score_df[base::rev(base::order(score_df$score_spearman)), ], ontology)
    enrichment_sobolev  <- UGenePred::enrichment_analysis(score_df[base::order(score_df$score_sobolev) , ], ontology)
    enrichment_fisher   <- UGenePred::enrichment_analysis(score_df[base::order(score_df$score_fisher) , ], ontology)

    # before binding, add column for showing what method result is originated from.
    enrichment_pearson[, "used_method"] <- base::rep("Pearson", base::nrow(enrichment_pearson))
    enrichment_spearman[, "used_method"] <- base::rep("Spearman", base::nrow(enrichment_spearman))
    enrichment_sobolev[, "used_method"] <- base::rep("Sobolev", base::nrow(enrichment_sobolev))
    enrichment_fisher[, "used_method"] <- base::rep("Fisher", base::nrow(enrichment_fisher))
    # bind all results by row
    combined_enrichment  <- base::rbind(enrichment_pearson, enrichment_spearman, enrichment_sobolev,  enrichment_fisher)

    # and keep the rows with the lowest corrected P-values.
    combined_enrichment  <- plyr::ddply(combined_enrichment, plyr::.(GOID), function(x) x[base::which.min(x$FDR),])
    # Order p-values on lowest > bigest
    enrichment_result   <- combined_enrichment[base::order(base::as.numeric(combined_enrichment$FDR)), ]
  }

  if(nrow(enrichment_result)>0) {
    # number the rownames and return the enrichment results.
    base::rownames(enrichment_result) <- c(1:base::nrow(enrichment_result))
    return(enrichment_result)
  } else {
    print("Could not find any similar genes!")
  }
}
