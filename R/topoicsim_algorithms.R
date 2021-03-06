#' GAPGOM internal - .topo_ic_sim_titj()
#'
#' This function is an internal function and should not be called by the user.
#'
#' Algorithm to calculate similarity between two GO terms.
#' This function is made for calculating topological similarity of two GO
#' terms in the GO DAG structure. The topological similarity is based on
#' edge weights and information content (IC). [1]
#'
#' @section Notes:
#' Internal function used in (topo_ic_sim_term()).
#'
#' @param go_id1 GO term of first term.
#' @param go_id2 GO term of second term.
#' @param topoargs list containing all the neccesary paramters/arguments for
#' the topoicsim algorithm, details can be view in the 
#' ".prepare_variables_topoicsim" function.
#' @return semantic topological similarity between the two given go terms. 
#'
#' @references [1] Ehsani R, Drablos F: \strong{TopoICSim: a new semantic
#' similarity measure based on gene ontology.} \emph{BMC Bioinformatics} 2016,
#' \strong{17}(1):296)
#'
#' @importFrom igraph igraph.from.graphNEL
#' @importFrom graph ftM2graphNEL subGraph
#' @importFrom AnnotationDbi get
#' @importFrom RBGL sp.between
#' @keywords internal
.topo_ic_sim_titj <- function(go_id1, go_id2, topoargs) {
  old <- options(stringsAsFactors = FALSE, warn = -1)
  on.exit(options(old), add = TRUE)
  # return 1 if the gos are the same.
  if (go_id1 == go_id2) {return(1)}
  common_ancestors <- .common_ancestors(go_id1, go_id2, topoargs$ontology, 
                                                        topoargs$organism,
                                        topoargs$go_annotation, topoargs$IC)
  D_ti_tj_x <- NULL
  if (length(common_ancestors) != 0 && !any(is.na(common_ancestors))) {
    # prepare subgraphs from two GO terms go_id1 and go_id2
    sg1 <- subGraph(c(get(go_id1, topoargs$go_annotation), go_id1),
                    topoargs$weighted_dag)
    sg1 <- .set_edge_weight(sg1, topoargs$IC)
    sg1 <- igraph.to.graphNEL(sg1)
    sg2 <- subGraph(c(get(go_id2, topoargs$go_annotation), go_id2),
                    topoargs$weighted_dag)
    sg2 <- .set_edge_weight(sg2, topoargs$IC)
    sg2 <- igraph.to.graphNEL(sg2)
    # also prepare goid ics
    IC_GOID_1 <- topoargs$IC[go_id1][[1]]
    IC_GOID_2 <- topoargs$IC[go_id2][[1]] 
    
    D_ti_tj_x <- lapply(common_ancestors, function(x) {
      # To identify all disjunctive common ancestors
      immediate_children_x <- switch(topoargs$ontology,
                                     MF = topoargs$GOMFCHILDREN[[x]],
                                     BP = topoargs$GOBPCHILDREN[[x]],
                                     CC = topoargs$GOCCCHILDREN[[x]])
      if (x != "all" & x != topoargs$root & !is.na(x) &
        length(intersect(immediate_children_x, common_ancestors)) == 0) {
        
        wLP_x_root <- .longest_path(topoargs$weighted_dag, 
          topoargs$go_annotation, x, topoargs$root, topoargs$IC)
        
        sp1 <- sp.between(sg1, go_id1, x)
        ic_sp1 <- sp1[[1]]$length
        length_sp1 <- length(sp1[[1]]$path_detail)
        sp2 <- sp.between(sg2, go_id2, x)
        ic_sp2 <- sp2[[1]]$length
        length_sp2 <- length(sp2[[1]]$path_detail)
        
        if (!is.na(IC_GOID_1) & IC_GOID_1 != 0)
          ic_sp1 <- ic_sp1 + (1/(2 * IC_GOID_1))
        if (!is.na(IC_GOID_2) & IC_GOID_2 != 0)
          ic_sp2 <- ic_sp2 + (1/(2 * IC_GOID_2))

        wSP_ti_tj_x <- (ic_sp1 + ic_sp2) * (length_sp1 + length_sp2 - 2)
        return(wSP_ti_tj_x/wLP_x_root)
      }
    })
    D_ti_tj_x <- unlist(D_ti_tj_x, FALSE, FALSE)
  }
  if (!is.null(D_ti_tj_x)) {
      sim <- round(1 - (atan(min(D_ti_tj_x, na.rm = TRUE))/(pi/2)), 3)
      return(sim)
  } else {
      return(0)
  }
}

#' GAPGOM internal - .topo_ic_sim_g1g2()
#'
#' This function is an internal function and should not be called by the user.
#'
#' Algorithm to calculate similarity between GO terms of two genes.
#' This function is made for calculating topological similarity of two Genes
#' given their GO terms in the GO DAG structure. The topological similarity is
#' based on edge weights and information content (IC). [1]
#'
#' @section Notes:
#' Internal function used in (topo_ic_sim_genes()).
#'
#' @param gene1 Gene ID of the first Gene.
#' @param gene2 Gene ID of the second Gene.
#' @param topoargs list containing all the neccesary paramters/arguments for
#' the topoicsim algorithm, details can be view in the 
#' ".prepare_variables_topoicsim" function.
#' @return List containing the following;
#' $GeneSim;
#' similarity between genes taken from the mean of all term 
#' similarities.
#' $AllGoPairs;
#' All possible GO combinations with their semantic distances (matrix)
#'
#' @references [1] Ehsani R, Drablos F: \strong{TopoICSim: a new semantic
#' similarity measure based on gene ontology.} \emph{BMC Bioinformatics} 2016,
#' \strong{17}(1):296)
#'
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @importFrom matrixStats rowMaxs colMaxs
#' @keywords internal
.topo_ic_sim_g1g2 <- function(gene1, gene2, topoargs) {
  old <- options(stringsAsFactors = FALSE, warn = -1)
  on.exit(options(old), add = TRUE)
  if (topoargs$debug) {
    message("\nWorking on genepair; ", gene1, ", ", gene2)
  }
  # get goids for both genes
  if (gene1 %in% names(topoargs$custom_genes1)) {
    gos1 <- topoargs$custom_genes1[[gene1]]
  } else if (gene1 %in% names(topoargs$custom_genes2)) {
    gos1 <- topoargs$custom_genes2[[gene1]]
  } else {
    gos1 <- as.character(topoargs$translation_to_goids[
      topoargs$translation_to_goids$ID==gene1,]$GO)
  }
  
  if (gene2 %in% names(topoargs$custom_genes1)) {
    gos2 <- topoargs$custom_genes1[[gene2]]
  } else if (gene1 %in% names(topoargs$custom_genes2)) {
    gos2 <- topoargs$custom_genes2[[gene2]]
  } else {
    gos2 <- as.character(topoargs$translation_to_goids[
      topoargs$translation_to_goids$ID==gene2,]$GO)
  }
  
  # return NA if goids sums are both 0 (no goids available to measure)
  if (sum(!is.na(gos1)) == 0 | sum(!is.na(gos2)) == 0) {
    return(list(GeneSim = NA, 
                AllGoPairs = topoargs$all_go_pairs))
  }
  # prepare score matrix
  scores <- .prepare_score_matrix_topoicsim(gos1, gos2)
  
  # get unique go pairs and loop through them.
  unique_pairs <- .unique_combos(gos1, gos2)
  for (i in seq_len(nrow(unique_pairs))) {
    pair <- unique_pairs[i,]
    go1 <- pair[[1]]
    go2 <- pair[[2]]
    if (topoargs$debug) {
      message("Working on gopair; ", go1, ", ", go2)
    }
    # grab gos out of all_gos (calculated before initialization)
    scores <- .set_values(go1, go2, scores, topoargs$all_go_pairs[go1, go2])
  }
  if (topoargs$garbage_collection) {
    gc()
  }
  # if score is NA, return.
  if (!sum(!is.na(scores))) {
    return(list(GeneSim = NA,
           AllGoPairs = topoargs$all_go_pairs))
  }
  # set all NA's to 0 for square root.
  scores[is.na(scores)] <- 0
  scores <- sqrt(scores)
  m <- length(gos1)
  n <- length(gos2)
  
  sim <- max(
    sum(rowMaxs(scores, rows = seq_len(m), na.rm = TRUE)
    )/m,
    sum(colMaxs(scores, cols = seq_len(n), na.rm = TRUE)
    )/n)
  
  sim <- round(sim, digits = 3)
  # return final score
  return(list(GeneSim = sim, 
              AllGoPairs = topoargs$all_go_pairs))
}

#' GAPGOM internal - .topo_ic_sim()
#'
#' This function is an internal function and should not be called by the user.
#'
#' Algorithm to calculate similarity between two gene vectors.
#' This function is made for calculating topological similarity between two
#' gene lists of which each gene has its GO terms in the GO DAG structure.
#' The topological similarity is based on edge weights and information
#' content (IC). The output it a nxn matrix depending on the vector lengths.
#' Intraset similarity can be calculated by comparing the same gene vector to
#' itself and using mean() on the output. The same can be done for Interset
#' similarity, but between two \strong{different} gene lists. [1]
#'
#' @section Notes:
#' Internal function used in (topo_ic_sim_genes()).
#'
#' @param gene_list1 The first gene vector of gene IDs. Note; THIS IS NOT THE
#' ENSEMBLID. Instead use the gene ID adopted by NCBI.
#' @param gene_list2 Same type as gene_list1, will be compared to gene_list1.
#' @param topoargs list containing all the neccesary paramters/arguments for
#' the topoicsim algorithm, details can be view in the 
#' ".prepare_variables_topoicsim" function.
#' 
#' @return List containing the following;
#' $GeneSim;
#' similarity between genes taken from the mean of all term 
#' similarities within those genes. (matrix of similarities between genes).
#' Take the mean of this matrix to either get the InterSetSim or IntraSetSim
#' depending on your input.
#' $AllGoPairs;
#' All possible GO combinations with their semantic distances (matrix)
#'  
#' @references [1] Ehsani R, Drablos F: \strong{TopoICSim: a new semantic
#' similarity measure based on gene ontology.} \emph{BMC Bioinformatics} 2016,
#' \strong{17}(1):296)
#' 
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @keywords internal
.topo_ic_sim_geneset <- function(gene_list1, gene_list2, topoargs) {
    old <- options(stringsAsFactors = FALSE, warn = -1)
    on.exit(options(old), add = TRUE)
    
    # set up score matrix in advance
    score_matrix <- .prepare_score_matrix_topoicsim(gene_list1, gene_list2)
    
    # only loop through necesary vectors (unique pairs)
    unique_pairs <- .unique_combos(gene_list1, gene_list2)
    
    if (topoargs$progress_bar) {
      message("Merging gene(set) results...")
      pb <- txtProgressBar(min = 0, max = nrow(unique_pairs), style = 3)
    } else if (topoargs$verbose) {
      message("Merging gene(set) results...")
    }
    
    for (i in seq_len(nrow(unique_pairs))) {
      if (topoargs$progress_bar) {
        setTxtProgressBar(pb, getTxtProgressBar(pb)+1)
      }
      pair <- unique_pairs[i,]
      gene1 <- pair[[1]]
      gene2 <- pair[[2]]
      genepair_result <- .topo_ic_sim_g1g2(gene1,
                                          gene2, 
                                          topoargs)
      score_matrix <- .set_values(gene1, gene2, score_matrix, 
                                  genepair_result$GeneSim)
      topoargs$all_go_pairs <- genepair_result$AllGoPairs
    }
    if (topoargs$progress_bar) {
      message("\nDone!")
    } else if (topoargs$verbose) {
      message("Done!")
    }
    return(list(GeneSim = score_matrix, 
                AllGoPairs = topoargs$all_go_pairs))
}

#' GAPGOM - topo_ic_sim_genes()
#'
#' Algorithm to calculate similarity between GO terms of two genes/genelists.
#'
#' This function is made for calculating topological similarity between two
#' gene vectors of which each gene has its GO terms in the GO DAG structure.
#' The topological similarity is based on edge weights and information
#' content (IC). The output it a nxn matrix depending on the vector lengths.
#' Intraset similarity can be calculated by comparing the same gene vector to
#' itself and using mean() on the output. The same can be done for Interset
#' similarity, but between two \strong{different} gene lists (IntraSet and 
#' InterSet similarities are only applicable to gene sets). [1]
#'
#' @param organism organism where to be scanned genes reside in, this option
#' is neccesary to select the correct GO DAG. Options are based on the org.db
#' bioconductor package;
#' http://www.bioconductor.org/packages/release/BiocViews.html#___OrgDb
#' Following options are available: "fly", "mouse", "rat", "yeast",
#' "zebrafish", "worm", "arabidopsis", "ecolik12", "bovine", "canine",
#' "anopheles", "ecsakai", "chicken", "chimp", "malaria", "rhesus", "pig",
#' "xenopus".
#' @param ontology desired ontology to use for similarity calculations.
#' One of three;
#' "BP" (Biological process),
#' "MF" (Molecular function) or
#' "CC" (Cellular Component).
#' @param genes1 Gene ID(s) of the first Gene (vector).
#' @param genes2 Gene ID(s) of the second Gene (vector).
#' @param custom_genes1 Custom genes added to the first list, needs to be a 
#' named list with the name being the arbitrary ID and the value being a vector
#' of GO terms.
#' @param custom_genes2 same as custom_genes1 but added to second gene list.
#' @param verbose set to true for more informative/elaborate output.
#' @param debug verbosity for debugging.
#' @param progress_bar Whether to show the progress of the calculation 
#' (default = FALSE)
#' @param garbage_collection whether to do R garbage collection. This is
#' useful for very large calculations/datasets, as it might decrease ram usage.
#' This option might however increase calculation time slightly.
#' @param use_precalculation wheter to use precalculated score matrix or not.
#' This speeds up calculation for the most frequent GO terms.
#' Only available for human, mouse with ids entrez/ensembl. Default is False
#' because this is the safest and most accurate option. Every update of org.Db
#' libraries makes this matrix outdated, so use at your own risk.
#' @param drop vector of evidences in go data structure you want to
#' skip (see set_go_data).
#' @param all_go_pairs dataframe of GO Term pairs with a column
#' representing similarity between the two. You can add the dataframe from
#' previous runs to improve performance (only works if the last result has
#' at least part of the genes of the current run). You can also use it for 
#' pre-calculation and getting the results back in a fast manner.
#' @param idtype id type of the genes you specified. default="ENTREZID". To see
#' other options, enter empty string.
#' @param go_data prepared go_data, from the set_go_data function. It is
#' practically the same as in GOSemSim, but with a slightly nicer interface.
#'
#' @return List containing the following;
#' $GeneSim;
#' similarity between genes taken from the mean of all term similarities 
#' (single gene). Or a nxn matrix of gene similarities. Intraset similarity can
#' be calculated by comparing the same gene vector to itself and using mean() 
#' on the output. The same can be done for Interset similarity, but between 
#' two \strong{different} gene vectors (gene vector). ;
#' $AllGoPairs;
#' All possible GO combinations with their semantic distances (matrix). NAs 
#' might be present in the matrix, these are GO pairs that didn't occur.
#' 
#' @examples
#' # single gene mode
#' result <- GAPGOM::topo_ic_sim_genes("human", "MF", "218", "501")
#' 
#' # genelist mode
#' list1 <- c("126133","221","218","216","8854","220","219","160428","224",
#' "222","8659","501","64577","223","217","4329","10840","7915","5832")
#' # ONLY A PART OF THE GENELIST IS USED BECAUSE OF R CHECK TIME CONTRAINTS
#' result <- GAPGOM::topo_ic_sim_genes("human", "MF", list1[1:2], 
#'                                                    list1[1:2])
#'
#' # with custom gene
#' custom <- list(cus1=c("GO:0016787", "GO:0042802", "GO:0005524"))
#' result <- GAPGOM::topo_ic_sim_genes("human", "MF", "218", "501", 
#'                                     custom_genes1 = custom)
#'
#' @references [1] Ehsani R, Drablos F: \strong{TopoICSim: a new semantic
#' similarity measure based on gene ontology.} \emph{BMC Bioinformatics} 2016,
#' \strong{17}(1):296)
#'
#' @export
topo_ic_sim_genes <- function(organism, ontology, genes1, genes2,
  custom_genes1 = NULL, custom_genes2 = NULL, verbose = FALSE, debug = FALSE,
  progress_bar = TRUE, garbage_collection = FALSE, use_precalculation = FALSE,
  drop = NULL, all_go_pairs = NULL, idtype = "ENTREZID", go_data = NULL) {
  old <- options(stringsAsFactors = FALSE, warn = -1)
  on.exit(options(old), add = TRUE)
  starttime <- Sys.time()
  
  if (!(.check_organism(organism) &
        .check_ontology(ontology) &
        .check_ifclass(genes1, "character", "genes1") &
        .check_ifclass(genes2, "character", "genes2") &
        .check_custom_gene(custom_genes1) &
        .check_custom_gene(custom_genes2) &
        .check_ifclass(verbose, "logical", "verbose", accept_null = FALSE) &
        .check_ifclass(debug, "logical", "debug", accept_null = FALSE) &
        .check_ifclass(progress_bar, "logical", "progress_bar", 
                       accept_null = FALSE) &
        .check_ifclass(garbage_collection, "logical", "garbage_collection", 
                       accept_null = FALSE) &
        .check_ifclass(use_precalculation, "logical", "use_precalculation", 
                       accept_null = FALSE) &
        .check_ifclass(drop, "character", "drop") &
        (.check_ifclass(all_go_pairs, "matrix", all_go_pairs) |
         .check_ifclass(all_go_pairs, "Matrix", all_go_pairs)) &
        .check_ifclass(go_data, "GOSemSimDATA", "go_data") &
        .check_idtype(idtype, organism))
  ) {
    stop("Error: one or more arguments are faulty!")
  }
  # if everything is ok, start preparing...
  topoargs <- .prepare_variables_topoicsim(organism, ontology, genes1, genes2, 
    custom_genes1, custom_genes2, drop, verbose, debug, progress_bar, 
    use_precalculation, garbage_collection, all_go_pairs, keytype = idtype,
    go_data = go_data)
  # append gene_lists with custom genes if neccesary
  if(!is.null(topoargs$custom_genes1)) {
    genes1 <- c(names(topoargs$custom_genes1), genes1)
  }
  if(!is.null(topoargs$custom_genes2)) {
    genes2 <- c(names(topoargs$custom_genes2), genes2)
  }
  
  # resolve unique genes and gos, also update all_go_pairs
  # only loop through necesary vectors (unique pairs)
  unique_pairs_genes <- .unique_combos(genes1, genes2)
  unique_pairs_gos <- .resolve_genes_unique_gos(unique_pairs_genes, topoargs)
  # last check before calculating...
  if (!.check_ifclass(unique_pairs_gos, "data.table", "unique_pairs_gos", 
                     accept_null = FALSE) |
      nrow(unique_pairs_gos) == 0) {
    message("NO GO PAIRS FOUND! (using correct keytype?)")
    return(list(GeneSim = NULL, AllGoPairs = NULL))
  }
  # update all_go_terms with .all_go_similarities, this function resolves only
  # necessary go's and skips one's that are 0 or not relevant.
  topoargs$all_go_pairs <- .all_go_similarities(unique_pairs_gos, topoargs, 
    drop=NULL)
  if (is.null(topoargs$all_go_pairs)) {
    return(list(GeneSim = NULL, AllGoPairs = NULL))
  }
  
  if ((length(genes1)) == 1 &
      (length(genes2)) == 1) {
    # single gene topo
    result <- .topo_ic_sim_g1g2(genes1, genes2, topoargs)
  } else {
    # multi gene topo
    result <- .topo_ic_sim_geneset(genes1, genes2, topoargs)
  }
  if (verbose) {
    message("Calculation time (in seconds):")
    message(difftime(Sys.time(), starttime, units = "secs"))
  }
  return(result)
}

#' GAPGOM - topo_ic_sim_term()
#' 
#' Algorithm to calculate similarity between two GO terms.
#' 
#' This function is made for calculating topological similarity of two GO
#' terms in the GO DAG structure. The topological similarity is based on
#' edge weights and information content (IC). [1]
#'
#' @param organism where to be scanned genes reside in, this option
#' is neccesary to select the correct GO DAG. Options are based on the org.db
#' bioconductor package;
#' http://www.bioconductor.org/packages/release/BiocViews.html#___OrgDb
#' Following options are available: "fly", "mouse", "rat", "yeast",
#' "zebrafish", "worm", "arabidopsis", "ecolik12", "bovine", "canine",
#' "anopheles", "ecsakai", "chicken", "chimp", "malaria", "rhesus", "pig",
#' "xenopus".
#' @param ontology desired ontology to use for similarity calculations.
#' One of three;
#' "BP" (Biological process),
#' "MF" (Molecular function) or
#' "CC" (Cellular Component).
#' @param go1 GO term of first term.
#' @param go2 GO term of second term.
#' @param go_data prepared go_data, from the set_go_data function. It is
#' practically the same as in GOSemSim, but with a slightly nicer interface.
#' 
#' @return TopoICSim score between the two terms.
#'
#' @examples 
#' result <- topo_ic_sim_term("human", "MF", "GO:0018478", "GO:0047105")
#'
#' @references [1] Ehsani R, Drablos F: \strong{TopoICSim: a new semantic
#' similarity measure based on gene ontology.} \emph{BMC Bioinformatics} 2016,
#' \strong{17}(1):296)
#'
#' @export
topo_ic_sim_term <- function(organism, ontology, go1, go2, go_data = NULL) {
  if (!(.check_organism(organism) &
        .check_ontology(ontology) &
        .check_ifclass(go1, "character", "go1", accept_null = FALSE) &
        .check_ifclass(go2, "character", "go2", accept_null = FALSE) &
        .check_ifclass(go_data, "GOSemSimDATA", "go_data"))
  ) {
    stop("Error: one or more arguments are faulty!")
  }
  topoargs <- .prepare_variables_topoicsim(organism, 
                                           ontology, 
                                           go1, 
                                           go2, 
                                           term_only = TRUE,
                                           go_data = go_data)
  return(.topo_ic_sim_titj(go1, go2, topoargs))
}
