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
#' 
#' # grab 1000 random ids to filter out between 100 and max length
#' sort_list <- rownames(GAPGOM::ft5_example_data@assayData$exprs)[
#'   sample(100:nrow(ft5_example_data@assayData$exprs), 1000)]
#' gid <- "chr10:101953059..101953074,-"
#' GAPGOM::expression_prediction_function(gid, 
#'                                        GAPGOM::ft5_example_data, 
#'                                        sort_list, 
#'                                        "mouse", 
#'                                        "BP")
#' @importFrom plyr ddply .
#' @import Biobase
#' @export
expression_prediction_function <- function(gene_id,
                                expression_set,
                                id_select_vector,
                                organism,
                                ontology,
                                enrichment_cutoff = 250,
                                method = "combine") {
  old <- options(stringsAsFactors = FALSE)
  on.exit(options(old), add = TRUE)
  # prepare the data with some special operations/vars that are needed later
  expression_data_sorted <- expression_set@assayData$exprs[id_select_vector,]
  expression_data_sorted_ids <- rownames(expression_data_sorted)

  # Target expression data where gene id matches
  target_expression_data <- expression_set@assayData$exprs[gene_id,]
  
  # Generate the translation df using gosemsim.
  print("Looking up GO terms...")
  id_translation_df <- .generate_translation_df(expression_set, organism, 
                                                ontology)
  
  # make args list for ambiguous functions
  args <- list(
    "gene_id" = gene_id,
    "expression_set" = expression_set,
    "expression_data_sorted" = expression_data_sorted,
    "expression_data_sorted_ids" = expression_data_sorted_ids,
    "id_select_vector" = id_select_vector,
    "id_translation_df" = id_translation_df,
    "target_expression_data" = target_expression_data,
    "organism" = organism,
    "ontology" = ontology,
    "enrichment_cutoff" = enrichment_cutoff)


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
  if (nrow(enrichment_result) > 0) {
    # number the rownames and return the enrichment results.
    rownames(enrichment_result) <- c(1:nrow(enrichment_result))
    return(enrichment_result)
  } else {
    print("Could not find any similar genes!")
  }
}

#' GAPGOM internal - Ambiguous/prediction functions
#'
#' These functions are ambiguous/standardized functions that help make the
#' enrichment analysis steps more streamlined. Most measures have their own
#' specific way of their scores being calculated. For this, they have their own
#' function clauses. translation lookup from entrez to goid is also contained
#' within these functions.
#'
#' @section Notes:
#' These functions are internal functions and should not be called by the user.
#'
#' @name ambiguous_functions
NULL

#' @rdname ambiguous_functions
.ambiguous_scorecalc <- compiler::cmpfun(function(args, applyfunc) {
  # apply the score calculation function
  score <- apply(args$expression_data_sorted,
                       1, applyfunc)
  # prepare and format a datafrane to return
  score_df <- .prepare_score_df(args$expression_data_sorted_ids, score, args$gene_id)
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
    args$expression_set,
    args$id_select_vector,
    args$id_translation_df,
    args$organism,
    args$ontology,
    enrichment_cutoff = args$enrichment_cutoff
  )
  return(enrichment_result)
})

#' @rdname ambiguous_functions
.ambiguous_method_origin <- compiler::cmpfun(function(enrichment_result, 
                                                     methodname) {
  # add used_method column
  enrichment_result[, "used_method"] <- rep(
    methodname, nrow(enrichment_result))
  return(enrichment_result)
})

#' @rdname ambiguous_functions
.predict_pearson <- compiler::cmpfun(function(args) {
  score_df <- .ambiguous_scorecalc(args, function(x) abs(cor(
    as.numeric(x), args$target_expression_data)))
  enrichment_result <- .ambiguous_enrichment(args,
                                            .ambiguous_score_rev_sort(score_df))
  enrichment_result <- .ambiguous_method_origin(enrichment_result, "pearson")
  return(enrichment_result)
})

#' @rdname ambiguous_functions
.predict_spearman <- compiler::cmpfun(function(args) {
  score_df <- .ambiguous_scorecalc(args, function(x) abs(cor(
    as.numeric(x), args$target_expression_data, method = "spearman")))
  enrichment_result <- .ambiguous_enrichment(args,
                                            .ambiguous_score_rev_sort(score_df))
  enrichment_result <- .ambiguous_method_origin(enrichment_result, "spearman")
  return(enrichment_result)
})

#' @rdname ambiguous_functions
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
.predict_combined <- compiler::cmpfun(function(args) {
  # Run enrichment for each method
  enrichment_pearson <- .predict_pearson(args)
  enrichment_spearman <- .predict_spearman(args)
  enrichment_kendall <- .predict_kendall(args)
  enrichment_sobolev <- .predict_sobolev(args)
  enrichment_fisher <- .predict_fisher(args)
  
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

#' @rdname ambiguous_functions
.generate_translation_df <- function(expression_set, organism, ontology) {
  entrezid_col <- .resolve_entrezid_col(expression_set)
  go_data <- .set_go_data(organism, ontology, computeIC = F)
  trans_df <- new.env()
  trans_df$rowtracker <- 0
  trans_df$passed_ids <- list()
  entrez_go_dfs <- sapply(expression_set@featureData@data[,entrezid_col], function(entrezrawid) {
    trans_df$rowtracker <- trans_df$rowtracker + 1
    #print(rowtracker)
    #print(rownames(expression_set@assayData$exprs)[rowtracker])
    entrez_id <- unlist(strsplit(entrezrawid, ":"))[2]
    # test if id has already occured earlier
    goids <- trans_df$passed_ids[[entrez_id]]
    if (is.null(goids)) {
      goids <- go_data@geneAnno[go_data@geneAnno==entrez_id,]$GO # fix double lookup  
      trans_df$passed_ids[[entrez_id]] <- c(goids)
      }
    #View(goids)
    if (length(goids) != 0){
      return(data.frame(ORIGID=rownames(expression_set@assayData$exprs)[trans_df$rowtracker], ENTREZID=entrezrawid, GO=goids))
    }
  })
  entrez_go_df <- unique(as.data.frame(data.table::rbindlist(entrez_go_dfs)))
  return(entrez_go_df)
}

#' resolve function for entrez
.resolve_entrezid_col <- function(expression_set) {
  colnames_vector <- colnames(expression_set@featureData@data)
  exp <- regexec(".*entrez.*", colnames_vector)
  regex_result <- unlist(regmatches(colnames_vector, exp))
  if (length(regex_result) < 1) {
    return(NULL)
  } else {
    return(regex_result)
  }
}
