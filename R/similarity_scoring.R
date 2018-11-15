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
#' @param ncluster Amount of cores you want to run the combined method on. 
#' Default=1. Does not work for other methods.
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
expression_semantic_scoring <- function(gene_id,
                                        expression_set,
                                        method = "combine",
                                        ncluster = 1) {
  old <- options(stringsAsFactors = FALSE, warn=-1)
  on.exit(options(old), add = TRUE)
  starttime <- Sys.time()
  
  # prepare the data with some special operations/vars that are needed later
  gene_id <- as.character(gene_id)
  # Target expression data where gene id matches
  target_expression_data <- expression_set@assayData$exprs[gene_id,]
  
  # make args list for ambiguous functions
  args <- list(
    "gene_id" = gene_id,
    "expression_data" = expression_set@assayData$exprs,
    "target_expression_data" = target_expression_data,
    "ncluster" = ncluster)
  scores <- switch(method,
         "pearson" =  .score_pearson(args),
         "spearman" = .score_spearman(args),
         "kendall" = .score_kendall(args),
         "fisher" = .score_fisher(args),
         "sobolev" = .score_sobolev(args),
         "combine" = .score_combined(args),
         function() {
           cat("Selected method; \"", method, "\" does not exist! Refer to the ",
               "help page for options.", sep = "")
           return(NULL)
         }
  )
  return(scores)
}

#' @rdname ambiguous_functions
#' @importFrom stats cor
.score_pearson <- compiler::cmpfun(function(args) {
  score_df <- .ambiguous_scorecalc(args, args$expression_data, function(x) abs(cor(
    as.numeric(x), args$target_expression_data)))
  score_df <- .ambiguous_method_origin(score_df, "pearson")
  return(score_df)
})

#' @rdname ambiguous_functions
#' @importFrom stats cor
.score_spearman <- compiler::cmpfun(function(args) {
  score_df <- .ambiguous_scorecalc(args, args$expression_data, function(x) abs(cor(
    as.numeric(x), args$target_expression_data, method = "spearman")))
  score_df <- .ambiguous_method_origin(score_df, "spearman")
  return(score_df)
})

#' @rdname ambiguous_functions
#' @importFrom stats cor
.score_kendall <- compiler::cmpfun(function(args) {
  score_df <- .ambiguous_scorecalc(args, args$expression_data, function(x) abs(cor(
    as.numeric(x), args$target_expression_data, method = "kendall")))
  score_df <- .ambiguous_method_origin(score_df, "kendall")
  return(score_df)
})

#' @rdname ambiguous_functions
.score_fisher <- compiler::cmpfun(function(args) {
  score_df <- .ambiguous_scorecalc(args, args$expression_data, function(x) fisher_metric(
    as.numeric(x), args$target_expression_data))
  score_df <- .ambiguous_method_origin(score_df, "fisher")
  return(score_df)
})

#' @rdname ambiguous_functions
.score_sobolev <- compiler::cmpfun(function(args) {
  score_df <- .ambiguous_scorecalc(args, args$expression_data, function(x) sobolev_metric(
    as.numeric(x), args$target_expression_data))
  score_df <- .ambiguous_method_origin(score_df, "sobolev")
  return(score_df)
})

#' @rdname ambiguous_functions
#' @importFrom parallel makeCluster stopCluster
#' @importFrom doParallel registerDoParallel
#' @import foreach
.score_combined <- compiler::cmpfun(function(args) {
  # Run enrichment for each method
  score_methods <- c(.score_pearson, .score_spearman, .score_kendall, 
                       .score_sobolev, .score_fisher)
  # cl <- makeCluster(args$ncluster)
  # registerDoParallel(cl)
  # scores <- foreach(method=score_methods, .combine = 'rbind') %dopar% {
  #   return(method(args))
  # }
  # stopCluster(cl)
  
  # apply for non-parallel testing purposes
  scores <<-
    lapply(score_methods, function(method) {
    return(method(args))
      })
  scores <- do.call("rbind", scores)
  
  return(scores)
})