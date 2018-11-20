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
#' @importFrom GO.db GOMFCHILDREN GOBPCHILDREN GOCCCHILDREN GOMFPARENTS
#' GOBPPARENTS GOCCPARENTS GOMFANCESTOR GOBPANCESTOR GOCCANCESTOR
#' @importFrom igraph igraph.from.graphNEL
#' @importFrom graph ftM2graphNEL subGraph
#' @importFrom AnnotationDbi get
#' @importFrom RBGL sp.between
#' @keywords internal
.topo_ic_sim_titj <- compiler::cmpfun(function(go_id1,
                                               go_id2,
                                               topoargs) {
  old <- options(stringsAsFactors = FALSE, warn = -1)
  on.exit(options(old), add = TRUE)
  # return 1 if the gos are the same.
  if (go_id1 == go_id2) {return(1)}
  common_ancestors <- .common_ancestors(go_id1, go_id2, topoargs$ontology, 
                                                        topoargs$organism,
                                        topoargs$go_annotation, topoargs$IC)
  D_ti_tj_x <- NULL
  if (length(common_ancestors) != 0 && !is.na(common_ancestors)) {
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
                                     MF = GOMFCHILDREN[[x]],
                                     BP = GOBPCHILDREN[[x]],
                                     CC = GOCCCHILDREN[[x]])
      if (x != "all" & x != topoargs$root & !is.na(x) &
        length(intersect(immediate_children_x, common_ancestors)) == 0) {
        # Subgraph from a disjunctive common ancestor to root
        sglca <- subGraph(c(get(x, topoargs$go_annotation), x), 
                          topoargs$weighted_dag)
        sglca <- .set_edge_weight(sglca, topoargs$IC)
        
        wLP_x_root <- .longest_path(sglca, x, topoargs$root, topoargs$IC)
        
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
    D_ti_tj_x <- unlist(D_ti_tj_x, F, F)
  }
  if (!is.null(D_ti_tj_x)) {
      sim <- round(1 - (atan(min(D_ti_tj_x, na.rm = TRUE))/(pi/2)), 3)
      return(sim)
  } else {
      return(0)
  }
})

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
#' @keywords internal
.topo_ic_sim_g1g2 <- compiler::cmpfun(function(gene1,
                                               gene2,
                                               topoargs) {
  old <- options(stringsAsFactors = FALSE, warn = -1)
  on.exit(options(old), add = TRUE)
  if (topoargs$verbose) {
    message(paste0("\nWorking on genepair; ", gene1, ", ", gene2))
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
  if (sum(!is.na(gos1)) == 0 || sum(!is.na(gos2)) == 0) {
    return(list(GeneSim = NA, 
                AllGoPairs = topoargs$all_go_pairs))
  }
  scores <- .prepare_score_matrix_topoicsim(gos1, gos2)
  
  unique_pairs <- .unique_combos(gos1, gos2)
  
  if (topoargs$progress_bar) {
    pb <- txtProgressBar(min = 0, max = nrow(unique_pairs), style = 3)
  } 
  for (i in seq_len(nrow(unique_pairs))) {
    pair <- unique_pairs[i,]
    if (topoargs$progress_bar) {
      setTxtProgressBar(pb, i)
    }
    go1 <- pair[[1]]
    go2 <- pair[[2]]
    if (topoargs$verbose) {
      message(paste0("Working on gopair; ", go1, ", ", go2))
    }
    # if this is not the case (row is not present), then run topo_ic_sim 
    # between 2 terms.
    if (!is.na(topoargs$all_go_pairs[go1, go2])) {
      # set already existing value.
      scores <- .set_values(go1, go2, scores, topoargs$all_go_pairs[go1, go2])
    } else if(topoargs$use_precalculation && 
              go1 %in% colnames(topoargs$selected_freq_go_pairs) && 
              go2 %in% colnames(topoargs$selected_freq_go_pairs)) {
      # set precalculated value.
      scores <- .set_values(go1, go2, scores, 
                            topoargs$selected_freq_go_pairs[go1, go2])
      topoargs$all_go_pairs <- .set_values(go1, go2, topoargs$all_go_pairs, 
                                           topoargs$selected_freq_go_pairs[
                                             go1, go2])
    } else {
      score <-
        .topo_ic_sim_titj(go1,
                          go2,
                          topoargs)
      
      if (topoargs$garbage_collection) {
        # garbage collection (topo term algorithm uses a lot of ram) only done
        # every 500th item/at end of loop
        if ((i %% 500) == 0) {
          gc()
        }
      }
      scores <- .set_values(go1, go2, scores, score)
      topoargs$all_go_pairs <- .set_values(go1, go2, topoargs$all_go_pairs, score)
    }
  }
  if (topoargs$progress_bar) {
    cat("\n")
  }
  if (topoargs$garbage_collection) {
    gc()
  }
  # if score is NA, return.
  if (!sum(!is.na(scores))) {
  return(list(GeneSim = NA,
              AllGoPairs = topoargs$all_go_pairs))
  }
  scores <- sqrt(scores)
  m <- length(gos1)
  n <- length(gos2)
  sim <- max(
    sum(sapply(seq_len(m), function(x) {
      max(scores[x, ], na.rm = TRUE)
  }))/m, 
    sum(sapply(seq_len(n), function(x) {
      max(scores[, x], na.rm = TRUE)
  }))/n)
  sim <- round(sim, digits = 3)
  # return final score
  return(list(GeneSim = sim, 
              AllGoPairs = topoargs$all_go_pairs))
})

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
.topo_ic_sim_geneset <- compiler::cmpfun(function(gene_list1,
                                         gene_list2,
                                         topoargs) {
    old <- options(stringsAsFactors = FALSE, warn = -1)
    on.exit(options(old), add = TRUE)
    
    # append gene_lists with custom genes if neccesary
    if(!is.null(topoargs$custom_genes1)) {
      gene_list1 <- c(names(topoargs$custom_genes1), gene_list1)
    }
    if(!is.null(topoargs$custom_genes2)) {
      gene_list2 <- c(names(topoargs$custom_genes2), gene_list2)
    }
    
    # set up score matrix in advance
    score_matrix <- .prepare_score_matrix_topoicsim(gene_list1, gene_list2)
    
    # only loop through necesary vectors (unique pairs)
    unique_pairs <- .unique_combos(gene_list1, gene_list2)
    
    if (topoargs$progress_bar) {
      pb <- txtProgressBar(min = 0, max = nrow(unique_pairs), style = 3)
    }
    
    # make a copy of topo arguments to turn off progressbar for genelevel
    topoargs_gen <- topoargs
    topoargs_gen$progress_bar <- F
    
    for (i in seq_len(nrow(unique_pairs))) {
      pair <- unique_pairs[i,]
      if (topoargs$progress_bar) {
        setTxtProgressBar(pb, i)
      }
      gene1 <- pair[[1]]
      gene2 <- pair[[2]]
      genepair_result <- .topo_ic_sim_g1g2(gene1,
                                          gene2, 
                                          topoargs_gen)
      score_matrix <- .set_values(gene1, gene2, score_matrix, 
                                  genepair_result$GeneSim)
      topoargs_gen$all_go_pairs <- genepair_result$AllGoPairs
    }
    topoargs$all_go_pairs <- topoargs_gen$all_go_pairs
    if (topoargs$progress_bar) {
      cat("\n")
    }
    return(list(GeneSim=score_matrix, 
                AllGoPairs = topoargs$all_go_pairs))
})

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
#' @param ontology desired ontology to use for similarity calculations.
#' One of three;
#' "BP" (Biological process),
#' "MF" (Molecular function) or
#' "CC" (Cellular Component).
#' @param organism organism where to be scanned genes reside in, this option
#' is neccesary to select the correct GO DAG. Options are based on the org.db
#' bioconductor package;
#' http://www.bioconductor.org/packages/release/BiocViews.html#___OrgDb
#' Following options are available: "fly", "mouse", "rat", "yeast",
#' "zebrafish", "worm", "arabidopsis", "ecolik12", "bovine", "canine",
#' "anopheles", "ecsakai", "chicken", "chimp", "malaria", "rhesus", "pig",
#' "xenopus".
#' @param genes1 Gene ID(s) of the first Gene (vector).
#' @param genes2 Gene ID(s) of the second Gene (vector).
#' @param custom_genes1 Custom genes added to the first list, needs to be a 
#' named list with the name being the arbitrary ID and the value being a vector
#' of GO terms.
#' @param custom_genes2 same as custom_genes1 but added to second gene list.
#' @param verbose set to true for more informative/elaborate output.
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
#' @param drop vector of GOID you want to exclude from the analysis.
#' @param all_go_pairs dataframe of GO Term pairs with a column
#' representing similarity between the two. You can add the dataframe from
#' previous runs to improve performance (only works if the last result has
#' at least part of the genes of the current run). You can also use it for 
#' pre-calculation and getting the results back in a fast manner.
#' @param idtype id type of the genes you specified. default="ENTREZID"
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
#' result <- GAPGOM::topo_ic_sim_genes("human", "MF", "218", "501", drop = NULL)
#' # genelist mode
#' list1 <- c("126133","221","218","216","8854","220","219","160428","224",
#' "222","8659","501","64577","223","217","4329","10840","7915")
#' result <- GAPGOM::topo_ic_sim_genes("human", "MF", list1, list1, drop = NULL)
#'
#' # with custom gene
#' custom <- list(cus1=c("GO:0016787", "GO:0042802", "GO:0005524"))
#' result <- GAPGOM::topo_ic_sim_genes("human", "MF", "218", "501", 
#'                                     custom_genes1 = custom, drop = NULL)
#'
#' @references [1] Ehsani R, Drablos F: \strong{TopoICSim: a new semantic
#' similarity measure based on gene ontology.} \emph{BMC Bioinformatics} 2016,
#' \strong{17}(1):296)
#'
#' @export
topo_ic_sim_genes <- compiler::cmpfun(function(organism,
                              ontology,
                              genes1, 
                              genes2,
                              custom_genes1 = NULL,
                              custom_genes2 = NULL,
                              verbose = FALSE,
                              progress_bar = TRUE,
                              garbage_collection = FALSE,
                              use_precalculation = FALSE, 
                              drop = NULL,
                              all_go_pairs = NULL,
                              idtype = "ENTREZID",
                              go_data = NULL
                              ) {
  old <- options(stringsAsFactors = FALSE, warn = -1)
  on.exit(options(old), add = TRUE)
  starttime <- Sys.time()
  .topo_ic_sim_argcheck_genes(ontology, organism, genes1, genes2)
  # if everything is ok, start preparing...
  topoargs <- .prepare_variables_topoicsim(organism, 
                                           ontology, 
                                           genes1, 
                                           genes2, 
                                           custom_genes1,
                                           custom_genes2,
                                           drop,
                                           verbose,
                                           progress_bar, 
                                           use_precalculation,
                                           garbage_collection, 
                                           all_go_pairs, 
                                           keytype = idtype,
                                           go_data = go_data)
  if ((length(genes1)+length(custom_genes1)) == 1 &&
      (length(genes2)+length(custom_genes2)) == 1) {
    # single gene topo
    if (verbose) {
      result <- .topo_ic_sim_g1g2(genes1, genes2, topoargs)
    } else {
      result <- suppressMessages(.topo_ic_sim_g1g2(genes1, genes2, topoargs))    
    }
  } else {
    # multi gene topo
    if (verbose) {
      result <- .topo_ic_sim_geneset(genes1, genes2, topoargs)
    } else {
      result <- suppressMessages(.topo_ic_sim_geneset(genes1, genes2, topoargs))    
    }
  }
  if (verbose) {message(Sys.time()-starttime)}
  return(result)
})

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
#' @importFrom GO.db GOMFCHILDREN GOBPCHILDREN GOCCCHILDREN GOMFPARENTS
#' GOBPPARENTS GOCCPARENTS GOMFANCESTOR GOBPANCESTOR GOCCANCESTOR
#' @importFrom igraph igraph.from.graphNEL
#' @importFrom graph ftM2graphNEL subGraph
#' @importFrom AnnotationDbi get
#' @importFrom RBGL sp.between
#' @export
topo_ic_sim_term <- compiler::cmpfun(function(organism,
                                              ontology,
                                              go1, 
                                              go2,
                                              go_data = NULL) {
  topoargs <- .prepare_variables_topoicsim(organism, 
                                           ontology, 
                                           go1, 
                                           go2, 
                                           term_only=T,
                                           go_data = go_data)
  return(.topo_ic_sim_titj(go1, go2, topoargs))
})