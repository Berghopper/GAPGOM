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
#' @param method default is 'combine'
#' @return The resulting matrix with prediction of similar GO terms.
#' @examples
#' expression_data <- read.table('/media/casper/USB_ccpeters/internship_thesis/papers/pred_lncRNA/lncRNA2function_data.txt', sep='\t', head=TRUE)
#' ensembl_id_to_goid <- read.table('/media/casper/USB_ccpeters/internship_thesis/papers/pred_lncRNA/EG2GO.txt', sep='\t',head=TRUE)
#' UGenePred::prediction_function('ENSG00000228630', 'BP')
#'
#' @export
prediction_function <- function(gene_id, ontology, expression_data = UGenePred::expression_data, ensembl_to_go_id_conversion_df = UGenePred::ensembl_id_to_go_id, start_expression_col = 4,
    end_expression_col = 22, enrichment_cutoff = 250, method = "combine") {
    # prepare the data with some special operations/vars that are needed later
    expression_data_pc <- expression_data[(expression_data$GeneType == "protein_coding"), ]
    ensembl_id_pc <- expression_data_pc$GeneID

    # Target expression data where gene id matches
    target_expression_data <- expression_data[(expression_data[, 1] == gene_id), ]
    target_expression_data <- base::as.numeric(target_expression_data[1, c(4:end_expression_col)])

    # make args list for functions
    args <- list("expression_data_pc" = expression_data_pc,
         "start_expression_col" = start_expression_col,
         "end_expression_col" = end_expression_col,
         "target_expression_data" = target_expression_data,
         "ensembl_id_pc" = ensembl_id_pc,
         "gene_id" = gene_id,
         "ontology" = ontology,
         "ensembl_to_go_id_conversion_df" = ensembl_to_go_id_conversion_df,
         "enrichment_cutoff" = enrichment_cutoff)


    # these functions calculate score between target expression of target gene vs the rest of the desired protein coding genes.  this is either a correlation metric or a
    # geometrical metric.  after score is calculated, these scores are enriched.
    enrichment_result <- NULL
    if (method == "Pearson") {
        enrichment_result <- predict_pearson(args)

    } else if (method == "Spearman") {
        enrichment_result <- predict_spearman(args)

    } else if (method == "FisherMetric") {
        enrichment_result <- predict_fisher(args)

    } else if (method == "SobolevMetric") {
        enrichment_result <- predict_sobolev(args)

    } else if (method == "combine") {
        # Run enrichment for each method
        enrichment_pearson <- predict_pearson(args)
        enrichment_spearman <- predict_spearman(args)
        enrichment_sobolev <- predict_sobolev(args)
        enrichment_fisher <- predict_fisher(args)

        # before binding, add column for showing what method result is originated from.
        enrichment_pearson[, "used_method"] <- base::rep(
          "Pearson", base::nrow(enrichment_pearson)
          )
        enrichment_spearman[, "used_method"] <- base::rep(
          "Spearman", base::nrow(enrichment_spearman)
          )
        enrichment_sobolev[, "used_method"] <- base::rep(
          "Sobolev", base::nrow(enrichment_sobolev)
          )
        enrichment_fisher[, "used_method"] <- base::rep(
          "Fisher", base::nrow(enrichment_fisher)
          )

        # bind all results by row
        combined_enrichment <- base::rbind(enrichment_pearson,
                                           enrichment_spearman,
                                           enrichment_sobolev,
                                           enrichment_fisher)
        # and keep the rows with the lowest corrected P-values.
        combined_enrichment <- plyr::ddply(combined_enrichment,
                                           plyr::.(GOID), function(x)
                                             x[base::which.min(x$FDR), ]
                                          )
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

ambiguous_scorecalc <- function(args, applyfunc) {
  score <- base::apply(args$expression_data_pc
                       [, c(args$start_expression_col:
                              args$end_expression_col)],
                       1, applyfunc)
  score_df <- prepare_score_df(args$ensembl_id_pc, score, args$gene_id)

  return(score_df)
}

ambiguous_score_rev_sort <- function(score_df) {
  return(score_df[base::rev(base::order(score_df[, 2])), ])
}

ambiguous_score_sort <- function(score_df) {
  return(score_df[base::order(score_df[, 2]), ])
}


ambiguous_enrichment <- function(args, ordered_score_df) {
  enrichment_result <- enrichment_analysis(
    ordered_score_df,
    args$ontology,
    ensembl_id_pc = args$ensembl_id_pc,
    ensembl_to_go_id_conversion_df = args$ensembl_to_go_id_conversion_df,
    enrichment_cutoff = args$enrichment_cutoff
  )
  return(enrichment_result)
}

predict_pearson <- function(args) {
  score_df <- ambiguous_scorecalc(args, function(x) base::abs(stats::cor(
    base::as.numeric(x), args$target_expression_data)))
  enrichment_result <- ambiguous_enrichment(args,
                                            ambiguous_score_rev_sort(score_df))
  return(enrichment_result)
}

predict_spearman <- function(args) {
  score_df <- ambiguous_scorecalc(args, function(x) base::abs(stats::cor(
    base::as.numeric(x), args$target_expression_data, method = "spearman")))
  enrichment_result <- ambiguous_enrichment(args,
                                            ambiguous_score_rev_sort(score_df))
  return(enrichment_result)
}

predict_fisher <- function(args) {
  score_df <- ambiguous_scorecalc(args, function(x) UGenePred::fisher_metric(
    as.numeric(x), args$target_expression_data))
  enrichment_result <- ambiguous_enrichment(args,
                                            ambiguous_score_sort(score_df))
  return(enrichment_result)
}

predict_sobolev <- function(args) {
  score_df <- ambiguous_scorecalc(args, function(x) UGenePred::sobolev_metric(
    as.numeric(x), args$target_expression_data))
  enrichment_result <- ambiguous_enrichment(args,
                                            ambiguous_score_sort(score_df))
  return(enrichment_result)
}
