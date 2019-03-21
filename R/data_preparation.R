### TOPOICSIM FUNCTIONS

#' GAPGOM internal - .prepare_variables_topoicsim()
#'
#' This function is an internal function and should not be called by the user.
#'
#' Prepares a list of intermediary arguments/parameters needed by topocisim.
#'
#' @section Notes:
#' Internal function used in (topo_ic_sim_genes).
#'
#' @param ontology desired ontology to use for similarity calculations.
#' One of three;
#' "BP" (Biological process),
#' "MF" (Molecular function) or
#' "CC" (Cellular Component).
#' @param organism where to be scanned genes reside in, this option
#' is neccesary to select the correct GO DAG. Options are based on the org.db
#' bioconductor package;
#' http://www.bioconductor.org/packages/release/BiocViews.html#___OrgDb
#' Following options are available: "fly", "mouse", "rat", "yeast",
#' "zebrafish", "worm", "arabidopsis", "ecolik12", "bovine", "canine",
#' "anopheles", "ecsakai", "chicken", "chimp", "malaria", "rhesus", "pig",
#' "xenopus".
#' @param gene_list1 The first gene vector of gene IDs. Note; THIS IS NOT THE
#' ENSEMBLID. Instead use the gene ID adopted by NCBI.
#' @param gene_list2 Same type as gene_list1, will be compared to gene_list1.
#' @param drop vector of evidences in go data structure you want to
#' skip (see set_go_data).
#' @param verbose set to true for more informative/elaborate output.
#' @param progress_bar Whether to show the progress of the calculation 
#' (default = FALSE)
#' @param garbage_collection whether to do R garbage collection. This is
#' useful for very large calculations/datasets, as it might decrease ram usage.
#' This option might however increase calculation time.
#' @param all_go_pairs dataframe of GO Term pairs with a column
#' representing similarity between the two. You can add the dataframe from
#' previous runs to improve performance (only works if the last result has
#' at least part of the genes of the current run). You can also use it for 
#' pre-calculation and getting the results back in a fast manner.
#' @param topoargs topoargs list that needs correction, default is a brand new
#' list
#' @param term_only specify if you only want arguments/params for TERM level.
#' @param keytype keytype used in querying of godata
#' @param go_data prepared go_data, from the set_go_data function. It is
#' practically the same as in GOSemSim, but with a slightly nicer interface.
#' 
#' @return list with all topoicsim arguments (can differ depending on 
#' algorithm level gene/geneset/term); organism, ontology, verbose, IC 
#' (Information Content from the go_data/go consortium), weighted_dag (DAG 
#' with weighted nodes), go_annotation (GO annotation of the correct root 
#' node), root (root node of the GO tree (MF, BP or CC) in 
#' annotationdbi-database string form), drop, progress_bar, garbage_collection, 
#' all_go_pairs, selected_freq_go_pairs (precalculated common GO_pairs that 
#' might increase performance), translation_to_goids (translation dataframe 
#' between ID and GOID.)
#' 
#' @importFrom GO.db GOMFCHILDREN GOBPCHILDREN GOCCCHILDREN GOMFPARENTS
#' GOBPPARENTS GOCCPARENTS GOMFANCESTOR GOBPANCESTOR GOCCANCESTOR
#' @importFrom AnnotationDbi toTable
#' @importFrom graph ftM2graphNEL
#' @importFrom utils sessionInfo
#' @keywords internal
.prepare_variables_topoicsim <- function(organism, ontology, gene_list1 = NULL, 
  gene_list2 = NULL, custom_genes1 = NULL, custom_genes2 = NULL, drop = NULL,
  verbose = FALSE, debug = FALSE, progress_bar = NULL, 
  use_precalculation = TRUE, garbage_collection = NULL, all_go_pairs = NULL, 
  topoargs = list(), term_only = FALSE, keytype = "ENTREZID", go_data = NULL) {
  # first get term arguments
  topoargs$organism <- organism
  topoargs$ontology <- ontology
  topoargs$verbose <- verbose
  topoargs$debug <- debug
  topoargs$custom_genes1 <- custom_genes1
  topoargs$custom_genes2 <- custom_genes2
  topoargs$use_precalculation <- use_precalculation
  # prepare go.db lookup environments
  topoargs$GOMFANCESTOR <- list2env(as.list(GOMFANCESTOR))
  topoargs$GOBPANCESTOR <- list2env(as.list(GOBPANCESTOR))
  topoargs$GOCCANCESTOR <- list2env(as.list(GOCCANCESTOR))
  topoargs$GOMFCHILDREN <- list2env(as.list(GOMFCHILDREN))
  topoargs$GOBPCHILDREN <- list2env(as.list(GOBPCHILDREN))
  topoargs$GOCCCHILDREN <- list2env(as.list(GOCCCHILDREN))
  
  
  if (verbose) {
    message("Preparing topoICSim data...")
    message("Preparing term data.")
  }
  # go_data --> IC
  if (is.null(topoargs$IC)) {
    if (verbose) {
      if (is.null(go_data)) {
        go_data <- set_go_data(organism = organism, ontology = ontology, 
                               keytype = keytype)
      }
    } else {
      if (is.null(go_data)) {
        go_data <- suppressMessages(set_go_data(organism = organism, 
                                                ontology = ontology, 
                                                keytype = keytype))  
      }
    }
    topoargs$IC <- go_data@IC
  }
  
  # xx_parents --> weighted dag
  if (is.null(topoargs$weighted_dag)) {
    xx_parents <- switch(ontology, MF = toTable(GOMFPARENTS),
                         BP = toTable(GOBPPARENTS), CC = toTable(GOCCPARENTS))
    topoargs$weighted_dag <- ftM2graphNEL(as.matrix(xx_parents[, 
                                                               c(seq_len(2))]))
  }
  # go_annotation
  if (is.null(topoargs$go_annotation)) {
    topoargs$go_annotation <- switch(ontology, 
                                     MF = topoargs$GOMFANCESTOR, 
                                     BP = topoargs$GOBPANCESTOR,
                                     CC = topoargs$GOCCANCESTOR)
  }
  # root
  if (is.null(topoargs$root)) {
    topoargs$root <- switch(ontology, MF = "GO:0003674", BP = "GO:0008150",
                            CC = "GO:0005575")
  }
  if (verbose) {
    message("Preparing gene/geneset data...")
  }
  # get gene arguments (if neccesary)
  if (!term_only) {
    topoargs$drop <- drop
    topoargs$progress_bar <- progress_bar
    topoargs$garbage_collection <- garbage_collection
    
    # precalculated values
    if (is.null(topoargs$selected_freq_go_pairs) &
        topoargs$use_precalculation) {
      if (keytype == "ENTREZID") {
        tmpkeytype <- "ENTREZ"
      } else {
        tmpkeytype <- keytype
      }
      topoargs$selected_freq_go_pairs <- freq_go_pairs[[paste0(tmpkeytype, "_", 
                                                               ontology, "_", 
                                                               organism)]]
      # check if precalculated matrix is up to date (Human)
      err_msg_both <- paste0("WARNING:\n",
                        "Precalculated matrix is out of date (\"%s\")! It", 
                        " will be updated in upcoming versions (please refer ",
                        "to the github issue page and report this error if it ",
                        "is not present yet; ",
                        "https://github.com/Berghopper/GAPGOM/issues). \n\nOr ",
                        "you may possibly have outdated packages, ",
                        "please check if your version of \"%s\" is lower than ",
                        "in the precalculated matrix: (%s). ",
                        "Turning option off...")
      
      if (organism == "human") {
        if (freq_go_pairs$
            sessioninfo$otherPkgs$org.Hs.eg.db$Version
            != .get_package_version("org.Hs.eg.db")) {
          topoargs$use_precalculation <- FALSE
          message(sprintf(err_msg_both, "org.Hs.eg.db", "org.Hs.eg.db", 
                          freq_go_pairs$
                            sessioninfo$otherPkgs$org.Hs.eg.db$Version))
        }
      } else {
        if (freq_go_pairs$
            sessioninfo$otherPkgs$org.Mm.eg.db$Version
            != .get_package_version("org.Mm.eg.db")) {
          topoargs$use_precalculation <- FALSE
          message(sprintf(err_msg_both, "org.Mm.eg.db", "org.Mm.eg.db", 
                          freq_go_pairs$
                            sessioninfo$otherPkgs$org.Mm.eg.db$Version))
        }  
      }
    }
    # translation_to_goids
    if (is.null(topoargs$translation_to_goids)) {
      if (is.null(go_data)) {
        go_data <- set_go_data(organism = organism, ontology = ontology, 
                               computeIC = FALSE, keytype = keytype)
      }
      if ((is.null(gene_list1) & is.null(custom_genes1)) |
          (is.null(gene_list2) & is.null(custom_genes2))) {
        topoargs$translation_to_goids <- NULL
      } else {
        topoargs$translation_to_goids <- .go_ids_lookup(unique(c(gene_list1, 
                                                                 gene_list2)), 
                                                        go_data, 
                                                        custom_genes = c(
                                                          custom_genes1, 
                                                          custom_genes2),
                                                        drop = drop) 
      }
    }
    if (is.null(all_go_pairs)) {
      if (is.null(topoargs$all_go_pairs)) {
        go_unique_list <- unique(c(topoargs$translation_to_goids$GO,
                                 unlist(custom_genes1, FALSE, FALSE),
                                 unlist(custom_genes2, FALSE, FALSE)))
        topoargs$all_go_pairs <- .prepare_score_matrix_topoicsim(go_unique_list,
                                                                 go_unique_list)
      }
    } else {
      go_unique_list <- unique(c(rownames(topoargs$all_go_pairs), 
                                 colnames(topoargs$all_go_pairs), 
                                 topoargs$translation_to_goids$GO,
                                 unlist(custom_genes1, FALSE, FALSE),
                                 unlist(custom_genes2, FALSE, FALSE)))
      topoargs$all_go_pairs <- .prepare_score_matrix_topoicsim(go_unique_list,
                                                               go_unique_list,
                                                               old_scores = 
                                                                 all_go_pairs)
    }
  }
  
  return(topoargs)
}


#' GAPGOM - set_go_data()
#'
#' Sets GO data like GOSemSim (this function purely makes choosing datasets a 
#' little easier and prints available keytypes if specified incorrectly.)
#'
#' @section Notes:
#' Internal function used in multiple functions of topoICSim.
#' 
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
#' @param computeIC whether to compute Information Content.
#' @param keytype keytype used in querying of godata
#' 
#' @return return godata as from GoSemSim
#' 
#' @examples 
#' # set go data for human, MF ontology.
#' go_data <- GAPGOM::set_go_data("human", "MF")
#' 
#' @import org.Hs.eg.db
#' @import org.Mm.eg.db
#' @importFrom GOSemSim godata
#' @importFrom AnnotationDbi keytypes
#' @export
set_go_data <- function(organism, ontology, computeIC = TRUE,
  keytype="ENTREZID") {
  species <- switch(organism, human = "org.Hs.eg.db",
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
                    stop(
                      print(paste0("Error, invalid organism; \"" , organism , 
                                   "\"!"))
                    ))
  # load correct library for GO data to check/show keytypes
  # eval(parse(text=paste0("library(\"",species,"\")")))
  if (!(keytype %in% keytypes(eval(parse(text=species))))) {
    stop("FATAL; SPECIFIED KEYTYPE; \"", keytype,
                "\" IS NOT AVAILABLE. AVAILABLE KEYTYPES;\n", 
         paste0(keytypes(eval(parse(text=species))), collapse = ", "))
  }
  return(godata(species, ont = ontology, computeIC = computeIC, keytype = 
                  keytype))
}

###TOPOICSIM FUNCTIONS

#' GAPGOM internal - .prepare_score_matrix_topoicsim()
#'
#' This function is an internal function and should not be called by the user.
#'
#' Prepared score dataframe for the semantic measures
#'
#' @section Notes:
#' Internal function used in multiple functions of topoICSim.
#'
#' @param vec1 vector1 with arbitrary lookup names
#' @param vec2 vector2 with arbitrary lookup names.
#' @param sparse whether to implement Matrix or matrix (default=F --> matrix)
#' @param old_scores old score matrix that has overlapping values with 
#' currently generated score matrix.
#' 
#' @return The score matrix with names and NA's.
#' @importFrom Matrix Matrix
#' @importFrom fastmatch %fin%
#' @keywords internal
.prepare_score_matrix_topoicsim <- function(vec1, vec2, sparse = FALSE, 
  old_scores = FALSE) {
  if (sparse) {
    score_matrix <- Matrix(nrow = length(vec1), 
         ncol = length(vec2), 
         dimnames = list(
           vec1, vec2))
  } else {
    score_matrix <- matrix(nrow = length(vec1), 
                           ncol = length(vec2), 
                           dimnames = list(
                             vec1, vec2))
    }
  score_matrix <- .set_identical_items(score_matrix)
  
  if (old_scores) {
    ### Add old scores (only intersecting)
    
    # first check which rownames/colnames intersect
    intersecting_rn <- rownames(old_scores)[rownames(old_scores) %fin%
                                              rownames(score_matrix)]
    intersecting_cn <- colnames(old_scores)[colnames(old_scores)
                                            %fin% colnames(score_matrix)]
    # also do this for the transposed old score matrix (to get all intersecting
    # pairs)
    intersecting_rn_t <- rownames(t(old_scores))[rownames(t(old_scores)) %fin%
                                                   rownames(score_matrix)]
    intersecting_cn_t <- colnames(t(old_scores))[colnames(t(old_scores)) %fin%
                                                   colnames(score_matrix)]
    # set values
    score_matrix[intersecting_rn, intersecting_cn] <- old_scores[
      intersecting_rn, intersecting_cn]
    score_matrix[intersecting_rn_t, intersecting_cn_t] <- t(old_scores)[
      intersecting_rn_t, intersecting_cn_t]
  }
  return(score_matrix)
}

#' GAPGOM internal - .set_identical_items()
#'
#' This function is an internal function and should not be called by the user.
#'
#' Sets identical items within a score matrix, the similarity will always be 1.
#'
#' @section Notes:
#' Internal function used in various topoclsim algorithms.
#'
#' @param score_matrix score matrix for topoclsim
#' 
#' @return Same matrix with correct sets set to 1.
#' 
#' @importFrom fastmatch fmatch
#' @keywords internal
.set_identical_items <- function(score_matrix) {
  # set all matching names of matrix to 1 (Same genes).
  matched_by_row <- fmatch(rownames(score_matrix), colnames(score_matrix))
  for (row in which(!is.na(matched_by_row))) {
    col <- matched_by_row[row]
    score_matrix[row, col] <- 1.0
  }
  return(score_matrix)
}

#' GAPGOM internal - .unique_combos()
#'
#' This function is an internal function and should not be called by the user.
#'
#' Generates unique combinations of two vectors
#'
#' @section Notes:
#' Internal function used in various topoclsim algorithms.
#'
#' @param v1 first vector to be combined with second.
#' @param v2 second vector to be combined with first.
#' 
#' @return datatable with unique combinations of the vectors excluding items
#' with same name and double combinations
#' 
#' @import data.table
#' @importFrom plyr .
#' @keywords internal
.unique_combos <- function(v1, v2) {
  V1 <- NULL
  V2 <- NULL
  intermediate <- unique(CJ(v1, v2)[V1 > V2, c("V1", "V2") := .(V2, V1)])
  return(intermediate[V1 != V2])
}

### LNCRNAPRED FUNCTIONS

#' GAPGOM internal - .prepare_score_df()
#'
#' This function is an internal function and should not be called by the user.
#'
#' Prepared score dataframe for the enrichment analysis.
#'
#' @section Notes:
#' Internal function used in prediction_function().
#'
#' @param original_ids rowname ID's for selected expression data.
#' @param score score vector
#' @param gene_id Gene ID originally selected in the prediction function. 
#' (gets filtered)
#' 
#' @return The score dataframe with ensmbl ID's
#' 
#' @importFrom stats na.omit
#' @keywords internal
.prepare_score_df <- function(original_ids, score, gene_id) {
  # score_df is the 'score' (correlation/metric) dataframe against target gene.
  # combine the dataframe if combine is selected, otherwise just add regular 
  # score.
  score_df <- data.frame(original_ids, score)
  score_df <- na.omit(score_df)
  # filter gene ID's (Select everything except the chosen gene).
  score_df <- score_df[(score_df[, 1] != gene_id), ]
  return(score_df)
}