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
.go_ids_lookup <- compiler::cmpfun(function(ids, go_data, drop=NULL) {
  go_gene_anno <- data.table(go_data@geneAnno)
  rm(go_data)
  go_gene_anno <- go_gene_anno[!go_gene_anno$EVIDENCE %in% drop, 1:2]
  
  go_gene_anno <- unique(go_gene_anno[go_gene_anno$ENTREZID %in% ids,])
  
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
  rm(passed_ids)
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
#' @param ontology desired ontology to use for prediction. One of three;
#' "BP" (Biological process), "MF" (Molecular function) or "CC"
#' (Cellular Component). Cellular Component is not included with the package's
#' standard data and will thus yield no results. 
#' 
#' @return return the translation dataframe containing conversion from general 
#' ids to entrez/ensembl ids and goids.
#' 
#' @import data.table
.generate_translation_df <- compiler::cmpfun(function(expression_set, 
                                                      organism, 
                                                      ontology) {
  entrezid_col <- .resolve_entrezid_col(expression_set) # add keys support
  go_data <- .set_go_data(organism, ontology, computeIC = F)
  go_gene_anno <- unique(data.table(go_data@geneAnno)[,1:2])
  rm(go_data)
  
  # convert entrez_ids and grab subset of godata (quicker)
  all_entrezzes <- lapply(expression_set@featureData@data[,entrezid_col], 
                          function(entrezrawid) {
                            entrez_split <- unlist(strsplit(entrezrawid, 
                                                            ",|:"), F, F)
                            entrez_ids <- entrez_split[seq(2, 
                                                           length(entrez_split), 2)]
                            return(entrez_ids)
                          })
  all_entrezzes <- unique(unlist(all_entrezzes, F, F))
  # grab correct go data
  go_gene_anno <- unique(go_gene_anno[go_gene_anno$ENTREZID %in% 
                                        all_entrezzes,])
  
  # keep track of row to properly bind main ID
  rowtracker <- 0
  passed_ids <- list()
  
  entrez_go_dfs <- lapply(expression_set@featureData@data[,entrezid_col], 
                          function(entrezrawid) {
                            rowtracker <<- rowtracker + 1
                            entrez_split <- unlist(strsplit(entrezrawid, ",|:"))
                            entrez_ids <- entrez_split[seq(2, 
                                                           length(entrez_split), 
                                                           2)]
                            
                            # test if id has already occured earlier
                            goids <- passed_ids[[entrezrawid]]
                            if (is.null(goids)) {
                              goids <- go_gene_anno[go_gene_anno$ENTREZID %in% 
                                                      entrez_ids,]$GO
                              non_duplicated_goids <- goids[!duplicated(goids)]
                              passed_ids[[entrezrawid]] <<- 
                                c(non_duplicated_goids)
                            }
                            # check if an output exists, if so return.
                            if (length(goids) != 0){
                              return(CJ(ORIGID=rownames(
                                expression_set@assayData$exprs)[rowtracker], 
                                        ENTREZID=entrezrawid, GO=goids))
                            }
                          })
  # bind the results, filter uniques and return.
  entrez_go_df <- unique(as.data.frame(data.table::rbindlist(entrez_go_dfs)))
  return(entrez_go_df)
})

#' GAPGOM internal - .resolve_entrez_col() 
#'
#' This function is an internal function and should not be called by the user.
#'
#' resolves columnname for entrez ids within an expressoin_set
#' 
#' @section Notes:
#' Internal function used in .generate_translation_df().
#'
#' @param expression_set ExpressionSet object --> see Biobase package.
#' 
#' @return column name of entrez id
.resolve_entrezid_col <- compiler::cmpfun(function(expression_set) {
  colnames_vector <- colnames(expression_set@featureData@data)
  exp <- regexec(".*entrez.*", colnames_vector)
  regex_result <- unlist(regmatches(colnames_vector, exp), F, F)
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
.ext_id_to_term_id <- compiler::cmpfun(function(data, list_top_genes) {
  # add 1 for each go that is in list top genes.
  dtdata <- as.data.table(data)
  quantified_only <- dtdata[dtdata$ORIGID %in% as.data.table(list_top_genes)[[1]], .N, by=GO]
  rm(dtdata)
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
.term_id_to_ext_id <- compiler::cmpfun(function(data) {
  dtdata <- as.data.table(data)
  return(as.data.frame(dtdata[, .N, by=GO]))
})

