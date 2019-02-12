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
.go_ids_lookup <- function(ids, go_data, custom_genes=NULL, drop=NULL) {
  go_gene_anno <- data.table(go_data@geneAnno)
  go_gene_anno <- go_gene_anno[!go_gene_anno$EVIDENCE %in% drop, c(seq_len(2))]
  
  go_gene_anno <- unique(go_gene_anno[go_gene_anno[[1]] %in% ids,])
  
  passed_ids <- list()
  go_df <- data.frame()
  
  for (id in ids) {
    # test if id has already occured earlier
    goids <- passed_ids[[id]]
    if (is.null(goids)) {
      goids <- unique(go_gene_anno[go_gene_anno[[1]]==as.character(id),]$GO)
      passed_ids[[id]] <- c(goids)
    }
    if (length(goids) != 0) {
      go_df <- rbind(go_df, data.frame(ID=id, GO=goids))
    }
  }
  if (!is.null(custom_genes)) {
    go_df <- rbind(go_df, 
                data.table::rbindlist(lapply(names(custom_genes), 
                       function(id, cus_genes) {
                        return(data.frame(ID=id, GO=cus_genes[[id]]))
                       }, custom_genes)
                  ))
  }
  go_df <- unique(as.data.frame(go_df))
  return(go_df)
}

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
.set_values <- function(item1, item2, the_matrix, value) {
  # set opposite pair to the same value if it exists
  if (item1 %in% rownames(the_matrix) && item2 %in% colnames(the_matrix)) {
    the_matrix[item1, item2] <- value
  }
  if (item2 %in% rownames(the_matrix) && item1 %in% colnames(the_matrix)) {
    the_matrix[item2, item1] <- value 
  }
  return(the_matrix)
}

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
#' @importFrom Biobase featureData pData assayData
#' @keywords internal
.generate_translation_df <- function(expression_set, organism, ontology, 
  keytype, verbose = FALSE, go_data = NULL) {
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
  go_gene_anno <- unique(data.table(go_data@geneAnno)[,c(seq_len(2))])
  
  # convert entrez_ids and grab subset of godata (quicker)
  all_keys <- lapply(pData(featureData(expression_set))[, keys_col], 
                     .entrezraw_to_entrez)
  all_keys <- unique(unlist(all_keys, FALSE, FALSE))
  # grab correct go data
  go_gene_anno <- unique(go_gene_anno[go_gene_anno[[1]] %in% 
                                        all_keys,])
  # keep track of row to properly bind main ID
  rowtracker <- 0
  passed_ids <- list()
  id_go_df <- data.frame()
  
  for (i in seq_len(nrow(pData(featureData(expression_set))))) {
    rawid <- as.character(pData(featureData(expression_set))[i, keys_col])
    ids <- .entrezraw_to_entrez(rawid)
    # test if id has already occurred earlier
    goids <- passed_ids[[rawid]]
    if (is.null(goids)) {
      goids <- go_gene_anno[go_gene_anno[[1]] %in% 
                              ids,]$GO
      non_duplicated_goids <- goids[!duplicated(goids)]
      passed_ids[[rawid]] <- 
        c(non_duplicated_goids)
    }
    # check if an output exists, if so return.
    if (length(goids) != 0){
      id_go_df <- rbind(CJ(ORIGID=rownames(
        assayData(expression_set)[["exprs"]])[i], 
        ID=rawid, GO=goids), id_go_df)
    }
  }
  # bind the results, filter uniques and return.
  id_go_df <- unique(as.data.frame(id_go_df))
  return(id_go_df)
}

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
#'
#' @importFrom Biobase featureData pData assayData
#' @keywords internal
.resolve_keys_col <- function(expression_set, keytype) {
  colnames_vector <- colnames(pData(featureData(expression_set)))
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
}

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
#' @importFrom GO.db GO
#' @keywords internal
.ext_id_to_term_id <- function(data, list_top_genes) {
  # add 1 for each go that is in list top genes.
  dtdata <- as.data.table(data)
  # GO from GO.db isn't actually used because the keyword is not evaluated.
  # this import is just to circumvent bioccheck
  quantified_only <- dtdata[dtdata$ORIGID %in% 
                              as.data.table(list_top_genes)[[1]], .N, by=GO]
  non_quantified_gos <- unique(data[!(data$GO %in% quantified_only$GO), ]$GO)
  return(as.data.frame(rbind(quantified_only, 
                             list(non_quantified_gos, 
                                  rep(0, length(non_quantified_gos))))))
}

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
#' @importFrom GO.db GO
#' @keywords internal
.term_id_to_ext_id <- function(data) {
  dtdata <- as.data.table(data)
  # GO from GO.db isn't actually used because the keyword is not evaluated.
  # this import is just to circumvent bioccheck
  return(as.data.frame(dtdata[, .N, by=GO]))
}

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
.entrezraw_to_entrez <- function(rawid) {
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
}

.organism_to_species_lib <- function(organism) {
  return(switch(tolower(organism), human = "org.Hs.eg.db",
           fly = "org.Dm.eg.db",
           mouse = "org.Mm.eg.db",
           rat = "org.Rn.eg.db",
           yeast = "org.Sc.sgd.db",
           zebrafish = "org.Dr.eg.db",
           worm = "org.Ce.eg.db",
           arabidopsis = "org.At.tair.db",
           ecolik12 = "org.EcK12.eg.db",
           bovine = "org.Bt.eg.db",
           canine = "org.Cf.eg.db",
           anopheles = "org.Ag.eg.db",
           ecsakai = "org.EcSakai.eg.db",
           chicken = "org.Gg.eg.db",
           chimp = "org.Pt.eg.db",
           malaria = "org.Pf.plasmo.db",
           rhesus = "org.Mmu.eg.db",
           pig = "org.Ss.eg.db",
           xenopus = "org.Xl.eg.db",
           message("Error, invalid organism; \"", organism , "\"!")))
}

.get_package_version <- function(pckg_name) {
  inst_pckgs <- installed.packages()
  return(inst_pckgs[inst_pckgs[,"Package"] == pckg_name,"Version"])
}