###TOPOICSIM FUNCTIONS

#' GAPGOM internal - .go_ids_lookup()
#'
#' This function is an internal function and should not be called by the user.
#'
#' Looks up goids per id (entrez/ensembl).
#'
#' @section Notes:
#' Internal function used in ().
#'
#' @param ids general ids that you want to search for in godata.
#' @param go_data the queried godata neccesary for the lookup
#' @param drop list of evidences you want to ignore.
#' 
#' @return return the translation dataframe containing conversion from id to 
#' goids.
#' @import data.table
#' @keywords internal
.go_ids_lookup <- compiler::cmpfun(function(ids, go_data, drop=NULL) {
  go_gene_anno <- data.table(go_data@geneAnno)
  go_gene_anno <- go_gene_anno[!go_gene_anno$EVIDENCE %in% drop, 1:2]
  
  go_gene_anno <- unique(go_gene_anno[go_gene_anno[[1]] %in% ids,])
  
  passed_ids <- list()
  go_dfs <- lapply(ids, function(id) {
    # test if id has already occured earlier
    goids <- passed_ids[[id]]
    if (is.null(goids)) {
      goids <- unique(go_gene_anno[go_gene_anno[[1]]==as.character(id),]$GO)
      passed_ids[[id]] <<- c(goids)
    }
    if (length(goids) != 0) {
      return(data.frame(ID=id, GO=goids))
    }
  })
  go_df <- unique(data.table::rbindlist(go_dfs))
  return(go_df)
})

#' GAPGOM internal - .set_values()
#'
#' This function is an internal function and should not be called by the user.
#'
#' Sets values within a topoclsim score matrix.
#'
#' @section Notes:
#' Internal function used in expression_prediction_function().
#'
#' @param item1 first item
#' @param item2 second item
#' @param the_matrix score matrix
#' @param value value to be set
#' 
#' @return return the matrix with newly set items.
#' @keywords internal
.set_values <- compiler::cmpfun(function(item1, item2, the_matrix, value) {
  # set opposite pair to the same value if it exists
  if (item1 %in% rownames(the_matrix) && item2 %in% colnames(the_matrix)) {
    the_matrix[item1, item2] <- value
  } 
  if (item2 %in% rownames(the_matrix) && item1 %in% colnames(the_matrix)) {
    the_matrix[item2, item1] <- value 
  }
  return(the_matrix)
})

### LNCRNAPRED FUNCTIONS

#' GAPGOM internal - .generate_translation_df()
#'
#' This function is an internal function and should not be called by the user.
#'
#' Generates a translation df specific for ExpressionSets, this is then used
#' for looking up ids and their respective GOIDS.
#'
#' @section Notes:
#' Internal function used in expression_prediction_function().
#'
#' @param expression_set ExpressionSet object --> see Biobase package.
#' @param organism where to be scanned genes reside in, this option
#' is neccesary to select the correct GO DAG. Options are based on the org.db
#' bioconductor package;
#' http://www.bioconductor.org/packages/release/BiocViews.html#___OrgDb
#' Following options are available: "fly", "mouse", "rat", "yeast",
#' "zebrafish", "worm", "arabidopsis", "ecolik12", "bovine", "canine",
#' "anopheles", "ecsakai", "chicken", "chimp", "malaria", "rhesus", "pig",
#' "xenopus". Fantom5 data only has "human" and "mouse" available depending
#' on the dataset.
#' @param keytype keytype used in querying of godata/columnnames 
#' @param ontology desired ontology to use for prediction. One of three;
#' "BP" (Biological process), "MF" (Molecular function) or "CC"
#' (Cellular Component). Cellular Component is not included with the package's
#' standard data and will thus yield no results. 
#' @param verbose set to true for more informative/elaborate output.
#' 
#' @return return the translation dataframe containing conversion from general 
#' ids to entrez/ensembl ids (and others) and goids.
#' 
#' @import data.table
#' @keywords internal
.generate_translation_df <- compiler::cmpfun(function(expression_set, 
                                                      organism, 
                                                      ontology,
                                                      keytype,
                                                      verbose = FALSE,
                                                      go_data = NULL) {
  keys_col <- .resolve_keys_col(expression_set, keytype)
  if (is.null(go_data)) {
    if (verbose) {
      go_data <- set_go_data(organism, ontology, computeIC = FALSE, 
                             keytype = keytype)
    } else {
      go_data <- suppressMessages(set_go_data(organism, ontology, 
                                              computeIC = FALSE, 
                                              keytype = keytype))
    }  
  }
  go_gene_anno <- unique(data.table(go_data@geneAnno)[,1:2])
  
  # convert entrez_ids and grab subset of godata (quicker)
  all_keys <- lapply(expression_set@featureData@data[, keys_col], 
                     .entrezraw_to_entrez)
  all_keys <- unique(unlist(all_keys, FALSE, FALSE))
  # grab correct go data
  go_gene_anno <- unique(go_gene_anno[go_gene_anno[[1]] %in% 
                                        all_keys,])
  
  # keep track of row to properly bind main ID
  rowtracker <- 0
  passed_ids <- list()
  
  id_go_dfs <- lapply(expression_set@featureData@data[,keys_col], 
                          function(rawid) {
                            rowtracker <<- rowtracker + 1
                            rawid <- as.character(rawid)
                            ids <- .entrezraw_to_entrez(rawid)
                            # test if id has already occurred earlier
                            goids <- passed_ids[[rawid]]
                            if (is.null(goids)) {
                              goids <- go_gene_anno[go_gene_anno[[1]] %in% 
                                                      ids,]$GO
                              non_duplicated_goids <- goids[!duplicated(goids)]
                              passed_ids[[rawid]] <<- 
                                c(non_duplicated_goids)
                            }
                            # check if an output exists, if so return.
                            if (length(goids) != 0){
                              return(CJ(ORIGID=rownames(
                                expression_set@assayData$exprs)[rowtracker], 
                                        ID=rawid, GO=goids))
                            }
                          })
  # bind the results, filter uniques and return.
  id_go_df <- unique(as.data.frame(data.table::rbindlist(id_go_dfs)))
  return(id_go_df)
})

#' GAPGOM internal - .resolve_keys_col() 
#'
#' This function is an internal function and should not be called by the user.
#'
#' resolves columnname for arbitrary ids within an expression_set
#' 
#' @section Notes:
#' Internal function used in .generate_translation_df().
#'
#' @param expression_set ExpressionSet object --> see Biobase package.
#' @param keytype keytype used in querying of godata
#' 
#' @return column name of id
#' @keywords internal
.resolve_keys_col <- compiler::cmpfun(function(expression_set, keytype) {
  colnames_vector <- colnames(expression_set@featureData@data)
  if (keytype == "ENTREZID") {
    keytype <- "entrez"
  }
  # create regex
  regex_str <- vapply(strsplit(keytype, split="")[[1]], function(x) {
    return(paste0("[", toupper(x), tolower(x), "]"))}, character(1))
  regex_str <- paste0(".*", paste0(regex_str, collapse=""), ".*")
  exp <- regexec(regex_str, colnames_vector)
  # get result
  regex_result <- unlist(regmatches(colnames_vector, exp), FALSE, FALSE)
  if (length(regex_result) < 1) {
    return(NULL)
  } else {
    return(regex_result)
  }
})

#' GAPGOM internal - .ext_id_to_term_id()
#'
#' This function is an internal function and should not be called by the user.
#'
#' Quantifies amount of extracted go terms within the list select top genes.
#'
#' @section Notes:
#' Internal function used in enrichment_analysis().
#'
#' @param data extracted_genes with correct onto etc. (matrix)
#' @param list_top_genes The list_top_genes matrix
#' 
#' @return return the quantified matrix
#' 
#' @import data.table
#' @keywords internal
.ext_id_to_term_id <- compiler::cmpfun(function(data, list_top_genes) {
  # add 1 for each go that is in list top genes.
  dtdata <- as.data.table(data)
  quantified_only <- dtdata[dtdata$ORIGID %in% as.data.table(list_top_genes)[[1]], .N, by=GO]
  non_quantified_gos <- unique(data[!(data$GO %in% quantified_only$GO), ]$GO)
  return(as.data.frame(rbind(quantified_only, list(non_quantified_gos, rep(0, length(non_quantified_gos))))))
})

#' GAPGOM internal - .term_id_to_ext_id()
#'
#' This function is an internal function and should not be called by the user.
#'
#' Quantifies amount of total extracted go terms
#'
#' @section Notes:
#' Internal function used in enrichment_analysis().
#'
#' @param data extracted_genes with correct onto etc. (matrix)
#' 
#' @return return the quantified matrix
#' 
#' @import data.table
#' @keywords internal
.term_id_to_ext_id <- compiler::cmpfun(function(data) {
  dtdata <- as.data.table(data)
  return(as.data.frame(dtdata[, .N, by=GO]))
})

#' GAPGOM internal - .entrezraw_to_entrez()
#'
#' This function is an internal function and should not be called by the user.
#'
#' Converts a raw ID from fantom data to a normal ID.
#'
#' @section Notes:
#' Internal function used in enrichment_analysis().
#'
#' @param rawid Raw entrez id. E.G.; entrez:23456
#' 
#' @return return the normal ID
#' 
#' @keywords internal
.entrezraw_to_entrez <- compiler::cmpfun(function(rawid) {
  rawid <- as.character(rawid)
  
  # first check if the id is raw at all.
  regex_str <- "[Ee][Nn][Tt][Rr][Ee][Zz][Gg][Ee][Nn][Ee]:\\d*,?"
  exp <- regexec(regex_str, rawid)
  # get result
  regex_result <- unlist(regmatches(rawid, exp), FALSE, FALSE)
  # check if raw id
  if (length(regex_result) > 0) {
    # regex match!
    ids_split <- unlist(strsplit(rawid, 
                                 ",|:"), FALSE, FALSE)
    ids <- ids_split[seq(2, length(ids_split), 2)] 
  } else {
    # in any other case, or where the id is not entrez, return the input.
    ids <- rawid
  }
  return(ids)
})