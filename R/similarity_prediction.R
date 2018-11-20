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
#' @param go_amount minimal amount of gos that a result needs to have to be 
#' considered similar enough.
#' @param filter_pvals filters pvalues that are equal to 0 (Default=FALSE).
#' @param idtype idtype of the expression_data. If not correctly specified, 
#' error will specify available IDs. default="ENTREZID"
#' @param verbose set to true for more informative/elaborate output.
#' @param id_select_vector gene rowname(s) that you want to keep in the
#' dataset. For example, let's say you need to only include protein coding
#' genes. You then make a vector including only ids that are protein coding.
#' @param id_translation_df df with translations between ID and GOID. col1 = ID,
#' col2 = GOID.
#' @param go_data from set_go_data function. A GoSemSim go_data object.
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
#' # keep everything that is a protein coding gene
#' filter_vector <- GAPGOM::expset@featureData@data[(
#' GAPGOM::expset@featureData@data$GeneType=="protein_coding"),]$GeneID
#' # set gid and run.
#' gid <- "ENSG00000228630"
#' 
#' result <- GAPGOM::expression_prediction(gid, 
#'                                         GAPGOM::expset, 
#'                                         filter_vector, 
#'                                         "human", 
#'                                         "BP",
#'                                         id_translation_df = 
#'                                           GAPGOM::id_translation_df,
#'                                         method = "combine", verbose = TRUE, 
#'                                         filter_pvals = TRUE
#' )
#' @importFrom plyr ddply .
#' @import Biobase
#' @export
expression_prediction <- function(gene_id,
                                expression_set,
                                organism,
                                ontology,
                                enrichment_cutoff = 250,
                                method = "combine",
                                significance = 0.05,
                                go_amount = 5,
                                filter_pvals = FALSE,
                                idtype = "ENTREZID",
                                verbose = FALSE,
                                id_select_vector = NULL,
                                id_translation_df = NULL,
                                go_data = NULL) {
  old <- options(stringsAsFactors = FALSE, warn=2)
  on.exit(options(old), add = TRUE)
  starttime <- Sys.time()
  
  # prepare the data with some special operations/vars that are needed later
  
  gene_id <- as.character(gene_id)
  # check if selection vector is defined, if not, include everything.
  if (is.null(id_select_vector)) {
    id_select_vector <- rownames(expression_set@assayData$exprs)
  }
  
  id_select_vector
  expression_data_sorted <- expression_set@assayData$exprs[(rownames(
    expression_set@assayData$exprs) %in% id_select_vector),]
  expression_data_sorted_ids <- rownames(expression_data_sorted)
  
  # check if expression_data sorted has anything at all
  if (nrow(expression_data_sorted) == 0) {
    stop("Selection of id_select_vector is empty!")
  }

  # Target expression data where gene id matches
  target_expression_data <- expression_set@assayData$exprs[gene_id,]
  
  # Generate the translation df using gosemsim.
  if (verbose) {
    message("Looking up GO terms...")
  }
  
  if (is.null(id_translation_df)) {
    id_translation_df <- .generate_translation_df(
      expression_set, 
      organism, 
      ontology,
      idtype,
      verbose,
      go_data = go_data)
  }
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
    "go_amount" = go_amount,
    "filter_pvals" = filter_pvals)


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
    enrichment_result[factor_index] <- lapply(enrichment_result[factor_index], 
                                              as.character)
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
.ambiguous_score_rev_sort <- compiler::cmpfun(function(score_df) {
  # reverse sorts the score column
  if (length(score_df) == 0) {
    return(NULL)
  } else {
    return(score_df[rev(order(score_df[, 2])), ]) 
  }
})

#' @rdname ambiguous_functions
.ambiguous_score_sort <- compiler::cmpfun(function(score_df) {
  # sorts the score column
  if (length(score_df) == 0) {
    return(NULL)
  } else {
    return(score_df[order(score_df[, 2]), ])
  }
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
    args$filter_pvals,
    args$go_amount
  )
  return(enrichment_result)
})

#' @rdname ambiguous_functions
.ambiguous_method_origin <- compiler::cmpfun(function(df, 
                                                      methodname) {
  # add used_method column
  if (length(df) != 0) {
    df[, "used_method"] <- rep(
      methodname, nrow(df))
  }
  return(df)
})

#' @rdname ambiguous_functions
#' @importFrom stats cor
.predict_pearson <- compiler::cmpfun(function(args) {
  score_df <- .ambiguous_scorecalc(args, args$expression_data_sorted, 
                                   misc_pearson)
  enrichment_result <- .ambiguous_enrichment(args,
                                             .ambiguous_score_rev_sort(score_df))
  enrichment_result <- .ambiguous_method_origin(enrichment_result, "pearson")
  return(enrichment_result)
})

#' @rdname ambiguous_functions
#' @importFrom stats cor
.predict_spearman <- compiler::cmpfun(function(args) {
  score_df <- .ambiguous_scorecalc(args, args$expression_data_sorted, 
                                   misc_spearman)
  enrichment_result <- .ambiguous_enrichment(args,
                                             .ambiguous_score_rev_sort(score_df))
  enrichment_result <- .ambiguous_method_origin(enrichment_result, "spearman")
  return(enrichment_result)
})

#' @rdname ambiguous_functions
#' @importFrom stats cor
.predict_kendall <- compiler::cmpfun(function(args) {
  score_df <- .ambiguous_scorecalc(args, args$expression_data_sorted, 
                                   misc_kendall)
  enrichment_result <- .ambiguous_enrichment(args,
                                             .ambiguous_score_rev_sort(score_df))
  enrichment_result <- .ambiguous_method_origin(enrichment_result, "kendall")
  return(enrichment_result)
})

#' @rdname ambiguous_functions
.predict_fisher <- compiler::cmpfun(function(args) {
  score_df <- .ambiguous_scorecalc(args, args$expression_data_sorted, 
                                   misc_fisher)
  enrichment_result <- .ambiguous_enrichment(args,
                                             .ambiguous_score_sort(score_df))
  enrichment_result <- .ambiguous_method_origin(enrichment_result, "fisher")
  return(enrichment_result)
})

#' @rdname ambiguous_functions
.predict_sobolev <- compiler::cmpfun(function(args) {
  score_df <- .ambiguous_scorecalc(args, args$expression_data_sorted, 
                                   misc_sobolev)
  enrichment_result <- .ambiguous_enrichment(args,
                                             .ambiguous_score_sort(score_df))
  enrichment_result <- .ambiguous_method_origin(enrichment_result, "sobolev")
  return(enrichment_result)
})

#' @rdname ambiguous_functions
.predict_combined <- compiler::cmpfun(function(args) {
  # Run enrichment for each method
  predict_methods <- c(.predict_pearson, .predict_spearman, .predict_kendall, 
                       .predict_sobolev, .predict_fisher)
  combined_enrichment <-
    lapply(predict_methods, function(method) {
      return(method(args))
    })
  combined_enrichment <- do.call("rbind", combined_enrichment)

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


