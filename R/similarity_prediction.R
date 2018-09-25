#' GAPGOM - expression_prediction_function()
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
#' "expression_data" containing expression values. default is 4, catered to
#' the included dataset.
#' @param end_expression_col should be the last column inside "expression_data"
#'  containing expression values. default is 22, catered to the included
#'  dataset.
#' @param enrichment_cutoff cutoff number for the amount of genes to be
#' enriched in the enrichment analysis. (default is 250)
#' @param method which statistical method to use for the prediction, currently
#' there are 5 available; "pearson", "spearman", "kendall", "fisher", "sobolev"
#' and "combine".
#'
#' @return The resulting dataframe with prediction of similar GO terms.
#' These are ordered with respect to FDR values. The following columns will be
#' in the dataframe;
#' GOID - Gene Ontology ID,
#' Ontology - Ontology type (MF or BP),
#' FDR - False Positive Rate,
#' Term - description of GOID,
#' used_method - the used method to determine the ontology term similarity
#' @examples
#' # Example with default dataset, take a look at the data documentation
#' # to fully grasp what's going on.
#' GAPGOM::expression_prediction_function('ENSG00000228630', 'BP')
#'
#' @importFrom plyr ddply .
#' @export
expression_prediction_function <- compiler::cmpfun(function(gene_id,
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
       "expression_data" = expression_data,
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
  enrichment_result <- switch(method,
    "pearson" =  predict_pearson(args),
    "spearman" = predict_spearman(args),
    "kendall" = predict_kendall(args),
    "fisher" = predict_fisher(args),
    "sobolev" = predict_sobolev(args),
    "combine" = predict_combined(args),
    function() {
      cat("Selected method; \"", method, "\" does not exist! Refer to the ",
          "help page for options.", sep = "")
      return(NULL)
    }
  )
  if (nrow(enrichment_result) > 0) {
    # number the rownames and return the enrichment results.
    rownames(enrichment_result) <- c(1:nrow(enrichment_result))
    return(enrichment_result)
  } else {
    print("Could not find any similar genes!")
  }
})

#' GAPGOM internal - Ambiguous/prediction functions
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
ambiguous_scorecalc <- compiler::cmpfun(function(args, applyfunc) {
  # apply the score calculation function
  score <- apply(args$expression_data_pc
                       [, c(args$start_expression_col:
                              args$end_expression_col)],
                       1, applyfunc)
  # prepare and format a datafrane to return
  score_df <- prepare_score_df(args$ensembl_id_pc, score, args$gene_id)

  return(score_df)
})

#' @rdname ambiguous_functions
ambiguous_score_rev_sort <- compiler::cmpfun(function(score_df) {
  # reverse sorts the score column
  return(score_df[rev(order(score_df[, 2])), ])
})

#' @rdname ambiguous_functions
ambiguous_score_sort <- compiler::cmpfun(function(score_df) {
  # sorts the score column
  return(score_df[order(score_df[, 2]), ])
})

#' @rdname ambiguous_functions
ambiguous_enrichment <- compiler::cmpfun(function(args, ordered_score_df) {
  # run the enrichment analysis function.
  enrichment_result <- enrichment_analysis(
    ordered_score_df,
    args$ontology,
    expression_data = args$expression_data,
    ensembl_to_go_id_conversion_df = args$ensembl_to_go_id_conversion_df,
    enrichment_cutoff = args$enrichment_cutoff
  )
  return(enrichment_result)
})

#' @rdname ambiguous_functions
ambiguous_method_origin <- compiler::cmpfun(function(enrichment_result, 
                                                     methodname) {
  # add used_method column
  enrichment_result[, "used_method"] <- rep(
    methodname, nrow(enrichment_result))
  return(enrichment_result)
})

#' @rdname ambiguous_functions
predict_pearson <- compiler::cmpfun(function(args) {
  score_df <- ambiguous_scorecalc(args, function(x) abs(cor(
    as.numeric(x), args$target_expression_data)))
  enrichment_result <- ambiguous_enrichment(args,
                                            ambiguous_score_rev_sort(score_df))
  enrichment_result <- ambiguous_method_origin(enrichment_result, "pearson")
  return(enrichment_result)
})

#' @rdname ambiguous_functions
predict_spearman <- compiler::cmpfun(function(args) {
  score_df <- ambiguous_scorecalc(args, function(x) abs(cor(
    as.numeric(x), args$target_expression_data, method = "spearman")))
  enrichment_result <- ambiguous_enrichment(args,
                                            ambiguous_score_rev_sort(score_df))
  enrichment_result <- ambiguous_method_origin(enrichment_result, "spearman")
  return(enrichment_result)
})

#' @rdname ambiguous_functions
predict_kendall <- compiler::cmpfun(function(args) {
  score_df <- ambiguous_scorecalc(args, function(x) abs(cor(
    as.numeric(x), args$target_expression_data, method = "kendall")))
  enrichment_result <- ambiguous_enrichment(args,
                                            ambiguous_score_rev_sort(score_df))
  enrichment_result <- ambiguous_method_origin(enrichment_result, "kendall")
  return(enrichment_result)
})

#' @rdname ambiguous_functions
predict_fisher <- compiler::cmpfun(function(args) {
  score_df <- ambiguous_scorecalc(args, function(x) fisher_metric(
    as.numeric(x), args$target_expression_data))
  enrichment_result <- ambiguous_enrichment(args,
                                            ambiguous_score_sort(score_df))
  enrichment_result <- ambiguous_method_origin(enrichment_result, "fisher")
  return(enrichment_result)
})

#' @rdname ambiguous_functions
predict_sobolev <- compiler::cmpfun(function(args) {
  score_df <- ambiguous_scorecalc(args, function(x) sobolev_metric(
    as.numeric(x), args$target_expression_data))
  enrichment_result <- ambiguous_enrichment(args,
                                            ambiguous_score_sort(score_df))
  enrichment_result <- ambiguous_method_origin(enrichment_result, "sobolev")
  return(enrichment_result)
})

#' @rdname ambiguous_functions
predict_combined <- compiler::cmpfun(function(args) {
  # Run enrichment for each method
  enrichment_pearson <- predict_pearson(args)
  enrichment_spearman <- predict_spearman(args)
  enrichment_kendall <- predict_kendall(args)
  enrichment_sobolev <- predict_sobolev(args)
  enrichment_fisher <- predict_fisher(args)
  
  # bind all results by row
  combined_enrichment <- rbind(enrichment_pearson,
                               enrichment_spearman,
                               enrichment_kendall,
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
  return(enrichment_result)
})
