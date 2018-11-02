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
#' @param gene_id gene rowname to be compared to the other GO terms.
#' @param expression_set ExpressionSet class containing expression values and
#' other useful information, see GAPGOM::f5_example_data
#' documentation for further explanation of this type. If you want a custom
#' ExpressionSet you have to define one yourself.
#' @param id_select_vector gene rowname(s) that you want filtered out of the
#' dataset. For example, let's say you need to only include protein coding
#' genes. You then select all other genes that aren't and put the ids in this
#' vector.
#' @param organism where to be scanned genes reside in, this option
#' is neccesary to select the correct GO DAG. Options are based on the org.db
#' bioconductor package;
#' http://www.bioconductor.org/packages/release/BiocViews.html#___OrgDb
#' Following options are available: "fly", "mouse", "rat", "yeast",
#' "zebrafish", "worm", "arabidopsis", "ecolik12", "bovine", "canine",
#' "anopheles", "ecsakai", "chicken", "chimp", "malaria", "rhesus", "pig",
#' "xenopus". Fantom5 data only has "human" and "mouse" available depending
#' on the dataset.
#' @param ontology desired ontology to use for prediction. One of three;
#' "BP" (Biological process), "MF" (Molecular function) or "CC"
#' (Cellular Component). Cellular Component is not included with the package's
#' standard data and will thus yield no results.
#' @param enrichment_cutoff cutoff number for the amount of genes to be
#' enriched in the enrichment analysis. (default is 250)
#' @param method which statistical method to use for the prediction, currently
#' there are 5 available; "pearson", "spearman", "kendall", "fisher", "sobolev"
#' and "combine".
#' @param significance normalized p-values (fdr) that are below this number 
#' will be kept. has to be a float/double between 0-1. Default is 0.05
#' @param filter_pvals filters pvalues that are equal to 0 (Default=FALSE).
#' @param ncluster Amount of cores you want to run the combined method on. 
#' Default=1. Does not work for other methods.
#' @param idtype idtype of the expression_data. If not correctly specified, 
#' error will specify available IDs. default="ENTREZID"
#' @param verbose set to true for more informative/elaborate output.
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
#' # to fully grasp what's going on with making of the filter etc. (Biobase 
#' # ExpressionSet)
#' 
#' # make a filter of all the ids that contain a 0 in their expression row.
#' sort_list <- rownames(GAPGOM::ft5_example_data@assayData$exprs[
#' apply(GAPGOM::ft5_example_data@assayData$exprs, 1, 
#' function(row) {all(row==0)}),])
#' # set an artbitrary gene you want to find similarities for. (5th row in this
#' # case)
#' gid <- rownames(GAPGOM::ft5_example_data@assayData$exprs)[5]
#' result <- GAPGOM::expression_prediction_function(gid, 
#'                                        GAPGOM::ft5_example_data, 
#'                                        sort_list, 
#'                                        "mouse", 
#'                                        "BP"
#'                                        )
#' @importFrom plyr ddply .
#' @import Biobase
#' @export
expression_prediction_function <- function(gene_id,
                                expression_set,
                                id_select_vector,
                                organism,
                                ontology,
                                enrichment_cutoff = 250,
                                method = "combine",
                                significance = 0.05,
                                filter_pvals = FALSE,
                                ncluster = 1,
                                idtype = "ENTREZID",
                                verbose = F) {
  # prepare the data with some special operations/vars that are needed later
  old <- options(stringsAsFactors = FALSE, warn=-1)
  on.exit(options(old), add = TRUE)
  starttime <- Sys.time()
  
  expression_data_sorted <- expression_set@assayData$exprs[!(rownames(
    expression_set@assayData$exprs) %in% id_select_vector),]
  #if(expression_data_sorted)
  expression_data_sorted_ids <- rownames(expression_data_sorted)

  # Target expression data where gene id matches
  target_expression_data <- expression_set@assayData$exprs[gene_id,]
  
  # Generate the translation df using gosemsim.
  if (verbose) {
    message("Looking up GO terms...")
  }
  id_translation_df <- .generate_translation_df(
    expression_set, 
    organism, 
    ontology,
    idtype,
    verbose)
  
  # make args list for ambiguous functions
  args <- list(
    "gene_id" = gene_id,
    "expression_data_sorted" = expression_data_sorted,
    "expression_data_sorted_ids" = expression_data_sorted_ids,
    "id_select_vector" = id_select_vector,
    "id_translation_df" = id_translation_df,
    "target_expression_data" = target_expression_data,
    "organism" = organism,
    "ontology" = ontology,
    "enrichment_cutoff" = enrichment_cutoff,
    "significance" = significance,
    "filter_pvals" = filter_pvals,
    "ncluster" = ncluster)


  # these functions calculate score between target expression of target gene
  # vs the rest of the desired protein coding genes.  this is either a
  # correlation metric or a geometrical metric.
  # after score is calculated, these scores are enriched.
  enrichment_result <- NULL
  enrichment_result <- switch(method,
    "pearson" =  .predict_pearson(args),
    "spearman" = .predict_spearman(args),
    "kendall" = .predict_kendall(args),
    "fisher" = .predict_fisher(args),
    "sobolev" = .predict_sobolev(args),
    "combine" = .predict_combined(args),
    function() {
      cat("Selected method; \"", method, "\" does not exist! Refer to the ",
          "help page for options.", sep = "")
      return(NULL)
    }
  )
  if (length(enrichment_result) > 0 && nrow(enrichment_result) > 0) {
    # number the rownames and return the enrichment results.
    rownames(enrichment_result) <- c(1:nrow(enrichment_result))
    # set all factors to strings.
    factor_index <- sapply(enrichment_result, is.factor)
    enrichment_result[factor_index] <- lapply(enrichment_result[factor_index], as.character)
    if (verbose) {
      message("Calculation time (seconds):")
      message(Sys.time()-starttime)
    }
    return(enrichment_result)
  } else {
    message("Could not find any similar genes!")
  }
}

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
#' @keywords internal
NULL

#' @rdname ambiguous_functions
.ambiguous_scorecalc <- compiler::cmpfun(function(args, applyfunc) {
  # apply the score calculation function
  score <- apply(args$expression_data_sorted,
                       1, applyfunc)
  # prepare and format a datafrane to return
  score_df <- .prepare_score_df(args$expression_data_sorted_ids, score, 
                                args$gene_id)
  return(score_df)
})

#' @rdname ambiguous_functions
.ambiguous_score_rev_sort <- compiler::cmpfun(function(score_df) {
  # reverse sorts the score column
  return(score_df[rev(order(score_df[, 2])), ])
})

#' @rdname ambiguous_functions
.ambiguous_score_sort <- compiler::cmpfun(function(score_df) {
  # sorts the score column
  return(score_df[order(score_df[, 2]), ])
})

#' @rdname ambiguous_functions
.ambiguous_enrichment <- compiler::cmpfun(function(args, ordered_score_df) {
  # run the enrichment analysis function.
  enrichment_result <- .enrichment_analysis(
    ordered_score_df,
    args$id_select_vector,
    args$id_translation_df,
    args$organism,
    args$ontology,
    args$enrichment_cutoff,
    args$significance,
    args$filter_pvals
  )
  return(enrichment_result)
})

#' @rdname ambiguous_functions
.ambiguous_method_origin <- compiler::cmpfun(function(enrichment_result, 
                                                     methodname) {
  # add used_method column
  if (length(enrichment_result) != 0) {
    enrichment_result[, "used_method"] <- rep(
      methodname, nrow(enrichment_result))
  }
  return(enrichment_result)
})

#' @rdname ambiguous_functions
#' @importFrom stats cor
.predict_pearson <- compiler::cmpfun(function(args) {
  score_df <- .ambiguous_scorecalc(args, function(x) abs(cor(
    as.numeric(x), args$target_expression_data)))
  enrichment_result <- .ambiguous_enrichment(args,
                                            .ambiguous_score_rev_sort(score_df))
  enrichment_result <- .ambiguous_method_origin(enrichment_result, "pearson")
  return(enrichment_result)
})

#' @rdname ambiguous_functions
#' @importFrom stats cor
.predict_spearman <- compiler::cmpfun(function(args) {
  score_df <- .ambiguous_scorecalc(args, function(x) abs(cor(
    as.numeric(x), args$target_expression_data, method = "spearman")))
  enrichment_result <- .ambiguous_enrichment(args,
                                            .ambiguous_score_rev_sort(score_df))
  enrichment_result <- .ambiguous_method_origin(enrichment_result, "spearman")
  return(enrichment_result)
})

#' @rdname ambiguous_functions
#' @importFrom stats cor
.predict_kendall <- compiler::cmpfun(function(args) {
  score_df <- .ambiguous_scorecalc(args, function(x) abs(cor(
    as.numeric(x), args$target_expression_data, method = "kendall")))
  enrichment_result <- .ambiguous_enrichment(args,
                                            .ambiguous_score_rev_sort(score_df))
  enrichment_result <- .ambiguous_method_origin(enrichment_result, "kendall")
  return(enrichment_result)
})

#' @rdname ambiguous_functions
.predict_fisher <- compiler::cmpfun(function(args) {
  score_df <- .ambiguous_scorecalc(args, function(x) fisher_metric(
    as.numeric(x), args$target_expression_data))
  enrichment_result <- .ambiguous_enrichment(args,
                                            .ambiguous_score_sort(score_df))
  enrichment_result <- .ambiguous_method_origin(enrichment_result, "fisher")
  return(enrichment_result)
})

#' @rdname ambiguous_functions
.predict_sobolev <- compiler::cmpfun(function(args) {
  score_df <- .ambiguous_scorecalc(args, function(x) sobolev_metric(
    as.numeric(x), args$target_expression_data))
  enrichment_result <- .ambiguous_enrichment(args,
                                            .ambiguous_score_sort(score_df))
  enrichment_result <- .ambiguous_method_origin(enrichment_result, "sobolev")
  return(enrichment_result)
})

#' @rdname ambiguous_functions
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @import foreach
.predict_combined <- compiler::cmpfun(function(args) {
  # Run enrichment for each method
  predict_methods <- c(.predict_pearson, .predict_spearman, .predict_kendall, 
                       .predict_sobolev, .predict_fisher)
  cl <- makeCluster(args$ncluster)
  registerDoParallel(cl)
  combined_enrichment <- foreach(method=predict_methods, .combine = 'rbind'
                                 ) %dopar% {
    return(method(args))
  }
  stopCluster(cl)
  if (length(combined_enrichment) == 0) {
    return(combined_enrichment)
  }
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


