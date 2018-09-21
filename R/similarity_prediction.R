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

#' UGenePred - expression_prediction_function()
#'
#' Predicts annotation of un-annotated genes based on existing
#' Gene Ontology annotation data and correlated expression patterns.
#'
#' This function is specifically made for predicting lncRNA annotation by
#' assuming "guilt by association". For instance, the expression data in this
#' package is actually based on mRNA expression data, but correlated with
#' lncRNA. This expression data is the used in combination with mRNA GO
#' annotation to calculate similarity scores between GO terms,
#'
#' @param gene_id gene EnsemblID to be compared to the other GO terms.
#' @param ontology desired ontology to use for prediction. One of three;
#' "BP" (Biological process), "MF" (Molecular function) or "CC"
#' (Cellular Component). Cellular Component is not included with the package's
#' standard data and will thus yield no results.
#' @param expression_data expression data dataframe with at least 1 or more
#' rows. Some of the columns should be named certain ways and the dataframe
#' should contain certain information neccesary for calculation;
#' 1st neccesary col; name; "GeneID" (EnsemblIDs should go here).
#' 2nd neccesary col; name; "GeneType" (genetype as described in the Ensembl
#' database).
#' Next, a range of n-columns long, containing expression values, names do not
#' matter here. (You can use the existing dataset as an example)
#' @param ensembl_to_go_id_conversion_df dataframe with corresponding GOIDs and
#' EnsemblIDs. As well as the expression value dataframe, this dataframe should
#' contain certain information with specific columnnames as well;
#' 1st neccesary col; name; "ensembl_gene_id" (Ensembl gene IDs).
#' 2nd neccesary col; go_id; "go_id" (Gene Ontology term ID).
#' 3rd and last neccesary col; "namespace_1003" (Toplevel GO term).
#' (You can use the existing dataset as an example)
#' @param start_expression_col should be the starting column inside
#' "expression_data" containing expression values.
#' @param end_expression_col should be the last column inside "expression_data"
#'  containing expression values.
#' @param enrichment_cutoff default is 250
#' @param method default is 'combine'
#' @return The resulting dataframe with prediction of similar GO terms.
#' @examples
#' UGenePred::expression_prediction_function('ENSG00000228630', 'BP')
#'
#'
#' @importFrom plyr ddply .
#' @export
expression_prediction_function <- function(gene_id,
                                ontology,
                                expression_data =
                                  UGenePred::expression_data,
                                ensembl_to_go_id_conversion_df =
                                  UGenePred::ensembl_id_to_go_id,
                                start_expression_col = 4,
                                end_expression_col = 22,
                                enrichment_cutoff = 250,
                                method = "combine") {
    old <- options(stringsAsFactors = FALSE)
    on.exit(options(old), add = TRUE)
    # prepare the data with some special operations/vars that are needed later
    expression_data_pc <- expression_data[
      (expression_data$GeneType == "protein_coding"), ]
    ensembl_id_pc <- expression_data_pc$GeneID

    # Target expression data where gene id matches
    target_expression_data <- expression_data[
      (expression_data[, 1] == gene_id), ]
    target_expression_data <- as.numeric(
      target_expression_data[1, c(4:end_expression_col)])

    # make args list for ambiguous functions
    args <- list("expression_data_pc" = expression_data_pc,
         "start_expression_col" = start_expression_col,
         "end_expression_col" = end_expression_col,
         "target_expression_data" = target_expression_data,
         "ensembl_id_pc" = ensembl_id_pc,
         "gene_id" = gene_id,
         "ontology" = ontology,
         "ensembl_to_go_id_conversion_df" = ensembl_to_go_id_conversion_df,
         "enrichment_cutoff" = enrichment_cutoff)


    # these functions calculate score between target expression of target gene
    # vs the rest of the desired protein coding genes.  this is either a
    # correlation metric or a geometrical metric.
    # after score is calculated, these scores are enriched.
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

        # before binding, add column for showing what method result is
        # originated from.
        enrichment_pearson[, "used_method"] <- rep(
          "Pearson", nrow(enrichment_pearson)
          )
        enrichment_spearman[, "used_method"] <- rep(
          "Spearman", nrow(enrichment_spearman)
          )
        enrichment_sobolev[, "used_method"] <- rep(
          "Sobolev", nrow(enrichment_sobolev)
          )
        enrichment_fisher[, "used_method"] <- rep(
          "Fisher", nrow(enrichment_fisher)
          )

        # bind all results by row
        combined_enrichment <- rbind(enrichment_pearson,
                                           enrichment_spearman,
                                           enrichment_sobolev,
                                           enrichment_fisher)
        # and keep the rows with the lowest corrected P-values.
        combined_enrichment <- ddply(combined_enrichment,
                                           .(GOID), function(x)
                                             x[which.min(x$FDR), ]
                                          )
        # Order p-values on lowest > bigest
        enrichment_result <- combined_enrichment[order(
          as.numeric(combined_enrichment$FDR)), ]
    } else {
        cat("Selected method; \"", method, "\" does not exist! Refer to the ",
            "help page for options.", sep = "")
        return(NULL)
    }
    ###

    if (nrow(enrichment_result) > 0) {
        # number the rownames and return the enrichment results.
        rownames(enrichment_result) <- c(1:nrow(enrichment_result))
        return(enrichment_result)
    } else {
        print("Could not find any similar genes!")
    }
}

#' UGenePred internal - Ambiguous/prediction functions
#'
#' These functions are ambiguous/standardized functions that help make the
#' enrichment analysis steps more streamlined. Most measures have their own
#' specific way of their scores being calculated. For this, they have their own
#' function clauses.
#'
#' @section Notes:
#' These functions are internal functions and should not be called by the user.
#'
#' @name ambiguous_functions
NULL

#' @rdname ambiguous_functions
ambiguous_scorecalc <- function(args, applyfunc) {
  # apply the score calculation function
  score <- apply(args$expression_data_pc
                       [, c(args$start_expression_col:
                              args$end_expression_col)],
                       1, applyfunc)
  # prepare and format a datafrane to return
  score_df <- prepare_score_df(args$ensembl_id_pc, score, args$gene_id)

  return(score_df)
}

#' @rdname ambiguous_functions
ambiguous_score_rev_sort <- function(score_df) {
  # reverse sorts the score column
  return(score_df[rev(order(score_df[, 2])), ])
}

#' @rdname ambiguous_functions
ambiguous_score_sort <- function(score_df) {
  # sorts the score column
  return(score_df[order(score_df[, 2]), ])
}

#' @rdname ambiguous_functions
ambiguous_enrichment <- function(args, ordered_score_df) {
  # run the enrichment analysis function.
  enrichment_result <- enrichment_analysis(
    ordered_score_df,
    args$ontology,
    ensembl_id_pc = args$ensembl_id_pc,
    ensembl_to_go_id_conversion_df = args$ensembl_to_go_id_conversion_df,
    enrichment_cutoff = args$enrichment_cutoff
  )
  return(enrichment_result)
}

#' @rdname ambiguous_functions
predict_pearson <- function(args) {
  score_df <- ambiguous_scorecalc(args, function(x) abs(cor(
    as.numeric(x), args$target_expression_data)))
  enrichment_result <- ambiguous_enrichment(args,
                                            ambiguous_score_rev_sort(score_df))
  return(enrichment_result)
}

#' @rdname ambiguous_functions
predict_spearman <- function(args) {
  score_df <- ambiguous_scorecalc(args, function(x) abs(cor(
    as.numeric(x), args$target_expression_data, method = "spearman")))
  enrichment_result <- ambiguous_enrichment(args,
                                            ambiguous_score_rev_sort(score_df))
  return(enrichment_result)
}

#' @rdname ambiguous_functions
predict_fisher <- function(args) {
  score_df <- ambiguous_scorecalc(args, function(x) fisher_metric(
    as.numeric(x), args$target_expression_data))
  enrichment_result <- ambiguous_enrichment(args,
                                            ambiguous_score_sort(score_df))
  return(enrichment_result)
}

#' @rdname ambiguous_functions
predict_sobolev <- function(args) {
  score_df <- ambiguous_scorecalc(args, function(x) sobolev_metric(
    as.numeric(x), args$target_expression_data))
  enrichment_result <- ambiguous_enrichment(args,
                                            ambiguous_score_sort(score_df))
  return(enrichment_result)
}
