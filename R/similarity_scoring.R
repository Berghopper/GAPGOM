#' GAPGOM - expression_prediction()
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
#' # to fully grasp what's going on with making of the filter etc. (Biobase
#' # ExpressionSet)
#' 
#' # set an artbitrary gene you want to find similarities for. (5th row in this
#' # case)
#' gid <- "ENSG00000228630"
#' result <- GAPGOM::expression_semantic_scoring(gid,
#'                                               GAPGOM::expset)
#'
#' @import Biobase
#' @export
expression_semantic_scoring <- function(gene_id, expression_set,
  method = "combine") {
  old <- options(stringsAsFactors = FALSE, warn=-1)
  on.exit(options(old), add = TRUE)
  starttime <- Sys.time()
  
  # check inputs
  if (!(.check_ifclass(gene_id, "character", "gene_id", accept_null = FALSE) &&
        .check_ifclass(expression_set, "ExpressionSet", "expression_set",
                       accept_null = FALSE) &&
        .check_method(method)
  )) {
    stop("Error: one or more arguments are faulty!")
  }
  
  # prepare the data with some special operations/vars that are needed later
  gene_id <- as.character(gene_id)
  # Target expression data where gene id matches
  target_expression_data <- assayData(expression_set)[["exprs"]][gene_id,]
  
  # make args list for ambiguous functions
  args <- list(
    "gene_id" = gene_id,
    "expression_data" = assayData(expression_set)[["exprs"]],
    "target_expression_data" = target_expression_data)
  scores <- switch(method,
         "pearson" =  .score_pearson(args),
         "spearman" = .score_spearman(args),
         "kendall" = .score_kendall(args),
         "fisher" = .score_fisher(args),
         "sobolev" = .score_sobolev(args),
         "combine" = .score_combined(args),
         function() {
           cat("Selected method; \"", method, 
               "\" does not exist! Refer to the ",
               "help page for options.", sep = "")
           return(NULL)
         }
  )
  return(scores)
}

#' @rdname ambiguous_functions
.ambiguous_scorecalc <- function(args, expression_data, applyfunc) {
  # check if the expression_data has any 0 values in it and filter them.
  expression_data <- expression_data[!(apply(expression_data, 1,
                                          function(x) {all(x==0)})),]
  if (nrow(expression_data) == 0) {
   stop("dataset does not contain any non-zero expression values!")
  }
  # apply the score calculation function
  score <- apply(expression_data,
                 1, applyfunc, args)
  # prepare and format a dataframe to return
  score_df <- .prepare_score_df(rownames(expression_data), score, 
                                args$gene_id)
  return(score_df)
}

#' @rdname ambiguous_functions
#' @importFrom stats cor
.score_pearson <- function(args) {
  score_df <- .ambiguous_scorecalc(args, args$expression_data, misc_pearson)
  score_df <- .ambiguous_method_origin(score_df, "pearson")
  return(score_df)
}

#' @rdname ambiguous_functions
#' @importFrom stats cor
.score_spearman <- function(args) {
  score_df <- .ambiguous_scorecalc(args, args$expression_data, misc_spearman)
  score_df <- .ambiguous_method_origin(score_df, "spearman")
  return(score_df)
}

#' @rdname ambiguous_functions
#' @importFrom stats cor
.score_kendall <- function(args) {
  score_df <- .ambiguous_scorecalc(args, args$expression_data, misc_kendall)
  score_df <- .ambiguous_method_origin(score_df, "kendall")
  return(score_df)
}

#' @rdname ambiguous_functions
.score_fisher <- function(args) {
  score_df <- .ambiguous_scorecalc(args, args$expression_data, misc_fisher)
  score_df <- .ambiguous_method_origin(score_df, "fisher")
  return(score_df)
}

#' @rdname ambiguous_functions
.score_sobolev <- function(args) {
  score_df <- .ambiguous_scorecalc(args, args$expression_data, misc_sobolev)
  score_df <- .ambiguous_method_origin(score_df, "sobolev")
  return(score_df)
}

#' @rdname ambiguous_functions
.score_combined <- function(args) {
  # Run enrichment for each method
  score_methods <- c(.score_pearson, .score_spearman, .score_kendall, 
                       .score_sobolev, .score_fisher)
  # apply for non-parallel testing purposes
  scores <-
    lapply(score_methods, function(method) {
    return(method(args))
      })
  scores <- do.call("rbind", scores)
  
  return(scores)
}
