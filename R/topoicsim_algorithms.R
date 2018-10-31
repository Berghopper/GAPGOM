#' GAPGOM internal - topo_ic_sim_titj()
#'
#' This function is an internal function and should not be called by the user.
#'
#' Algorithm to calculate similarity between two GO terms.
#' This function is made for calculating topological similarity of two GO
#' terms in the GO DAG structure. The topological similarity is based on
#' edge weights and information content (IC). [1]
#'
#' @param go_id1 GO term of first term.
#' @param go_id2 GO term of second term.
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
#' @param weighted_dag DAG with weighted nodes.
#' @param go_annotation GO annotation of the correct root node.
#' @param root root node of the GO tree (MF, BP or CC) in databse string form.
#' @param IC Information content from the go_data.
#' @param all_go_pairs dataframe of GO Term pairs with a column
#' representing similarity between the two. it is recommended to always keep
#' this on default unless you know what you are doing.
#' @return list of score and all_go_pairs, the latter of which will be 
#' reused.
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
.topo_ic_sim_titj <- compiler::cmpfun(function(go_id1,
                                               go_id2,
                                               topoargs) {
  old <- options(stringsAsFactors = FALSE, warn = -1)
  on.exit(options(old), add = TRUE)
  common_ancestors <- .common_ancestors(go_id1, go_id2, topoargs$ontology, 
                                                        topoargs$organism,
                                        topoargs$go_annotation, topoargs$IC)
  D_ti_tj_x <- NULL
  if (length(common_ancestors) != 0 && !is.na(common_ancestors)) {
    D_ti_tj_x <- lapply(common_ancestors, function(x) {
      # To identify all disjunctive common ancestors
      immediate_children_x <- switch(topoargs$ontology,
                                     MF = GOMFCHILDREN[[x]],
                                     BP = GOBPCHILDREN[[x]],
                                     CC = GOCCCHILDREN[[x]])
      
      if (x != "all" & x != topoargs$root & !is.na(x) &
        length(intersect(immediate_children_x, common_ancestors)) == 0) {
        # Subgraph from two GO terms go_id1 and go_id2
        sg1 <- subGraph(c(get(go_id1, topoargs$go_annotation), go_id1),
                        topoargs$weighted_dag)
        sg2 <- subGraph(c(get(go_id2, topoargs$go_annotation), go_id2),
                        topoargs$weighted_dag)
        # Subgraph from a disjunctive common ancestor to root
        sglca <- subGraph(c(get(x, topoargs$go_annotation), x), topoargs$weighted_dag)
        sglca <- .set_edge_weight(sglca, topoargs$IC)
        wLP_x_root <- .longest_path(sglca, x, topoargs$root, topoargs$IC)
        sg1 <- .set_edge_weight(sg1, topoargs$IC)
        sg1 <- igraph.to.graphNEL(sg1)
        sp1 <- sp.between(sg1, go_id1, x)
        ic_sp1 <- sp1[[1]]$length
        length_sp1 <- length(sp1[[1]]$path_detail)
        sg2 <- .set_edge_weight(sg2, topoargs$IC)
        sg2 <- igraph.to.graphNEL(sg2)
        sp2 <- sp.between(sg2, go_id2, x)
        ic_sp2 <- sp2[[1]]$length
        length_sp2 <- length(sp2[[1]]$path_detail)
        IC_GOID_1 <- topoargs$IC[go_id1][[1]]
        IC_GOID_2 <- topoargs$IC[go_id2][[1]]
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

#' GAPGOM - topo_ic_sim_g1g2()
#'
#' Algorithm to calculate similarity between GO terms of two genes.
#'
#' This function is made for calculating topological similarity of two Genes
#' given their GO terms in the GO DAG structure. The topological similarity is
#' based on edge weights and information content (IC). [1]
#'
#' @param gene1 Gene ID of the first Gene.
#' @param gene2 Gene ID of the second Gene.
#' @param progress_bar Whether to show the progress of the calculation 
#' (default = FALSE)
#' @param garbage_collection whether to do R garbage collection. This is
#' useful for very large calculations/datasets, as it might decrease ram usage.
#' This option might however increase calculation time.
#' @param go_data correct GO data from specefic species.
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
#' @param drop vector of GOID you want to exclude from the analysis.
#' @param translation_to_goids translational dataframe between ID and GOID.
#' @param all_go_pairs dataframe of GO Term pairs with a column
#' representing similarity between the two. it is recommended to always keep
#' this on default unless you know what you are doing.
#'
#' @return List containing the following;
#' $GeneSim;
#' similarity between genes taken from the mean of all term 
#' similarities.
#' $GO1;
#' vector of GO's from first gene
#' $GO2;
#' vector of GO's from second gene
#' $AllGoPairs;
#' All possible GO combinations with their semantic distances (matrix)
#' 
#' @examples
#' result <- GAPGOM::topo_ic_sim_g1g2("218", "501", ont="MF", organism="human", 
#'                                    drop=NULL)
#'
#' @references [1] Ehsani R, Drablos F: \strong{TopoICSim: a new semantic
#' similarity measure based on gene ontology.} \emph{BMC Bioinformatics} 2016,
#' \strong{17}(1):296)
#'
#' @importFrom AnnotationDbi toTable
#' @importFrom graph ftM2graphNEL
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
.topo_ic_sim_g1g2 <- compiler::cmpfun(function(gene1,
                                               gene2,
                                               topoargs) {
  old <- options(stringsAsFactors = FALSE, warn = -1)
  on.exit(options(old), add = TRUE)
  
  # get goids for both genes
  gos1 <- as.character(topoargs$translation_to_goids[topoargs$translation_to_goids$ID==gene1,]$GO)
  gos2 <- as.character(topoargs$translation_to_goids[topoargs$translation_to_goids$ID==gene2,]$GO)

  # return NA if goids sums are both 0 (no goids available to measure)
  if (sum(!is.na(gos1)) == 0 || sum(!is.na(gos2)) == 0) {
    return(list(GeneSim = NA, GO1 = gos1, GO2 = gos2, 
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
    # if this is not the case (row is not present), then run topo_ic_sim 
    # between 2 terms.
    #print(dput(topoargs$all_go_pairs))
    if (!is.na(topoargs$all_go_pairs[go1, go2])) {
      # set already existing value.
      scores <- .set_values(go1, go2, scores, topoargs$all_go_pairs[go1, go2])
    } else if(go1 %in% colnames(topoargs$selected_freq_go_pairs) && 
              go2 %in% colnames(topoargs$selected_freq_go_pairs)) {
      # set precalculated value.
      scores <- .set_values(go1, go2, scores, topoargs$selected_freq_go_pairs[go1, go2])
      topoargs$all_go_pairs <- .set_values(go1, go2, topoargs$all_go_pairs, topoargs$selected_freq_go_pairs[go1, go2])
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
  return(list(GeneSim = NA, GO1 = gos1, GO2 = gos2, AllGoPairs = 
                topoargs$all_go_pairs))
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
  return(list(GeneSim = sim, GO1 = gos1, GO2 = gos2, AllGoPairs = 
                topoargs$all_go_pairs))
})

#' GAPGOM - topo_ic_sim()
#'
#' Algorithm to calculate similarity between two gene vectors.
#'
#' This function is made for calculating topological similarity between two
#' gene lists of which each gene has its GO terms in the GO DAG structure.
#' The topological similarity is based on edge weights and information
#' content (IC). The output it a nxn matrix depending on the vector lengths.
#' Intraset similarity can be calculated by comparing the same gene vector to
#' itself and using mean() on the output. The same can be done for Interset
#' similarity, but between two \strong{different} gene lists. [1]
#'
#' @param gene_list1 The first gene vector of gene IDs. Note; THIS IS NOT THE
#' ENSEMBLID. Instead use the gene ID adopted by NCBI.
#' @param gene_list2 Same type as gene_list1, will be compared to gene_list1.
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
#' @param drop vector of GOID you want to exclude from the analysis.
#' @param progress_bar Whether to show the progress of the calculation 
#' (default = TRUE)
#' @param garbage_collection whether to do R garbage collection. This is
#' useful for very large calculations/datasets, as it might decrease ram usage.
#' This option might however increase calculation time.
#' 
#' @return List containing the following;
#' $GeneSim;
#' similarity between genes taken from the mean of all term 
#' similarities within those genes. (matrix of similarities between genes).
#' Take the mean of this matrix to either get the InterSetSim or IntraSetSim
#' depending on your input.
#' $GO1;
#' vector of GO's from first gene
#' $GO2;
#' vector of GO's from second gene
#' $AllGoPairs;
#' All possible GO combinations with their semantic distances (matrix)
#' @examples
#' list1 <- c("126133","221","218","216","8854","220","219","160428","224",
#' "222","8659","501","64577","223","217","4329","10840","7915")
#' result <- GAPGOM::topo_ic_sim(list1, list1, ont="MF", organism="human", 
#'                               drop=NULL)
#' 
#' @references [1] Ehsani R, Drablos F: \strong{TopoICSim: a new semantic
#' similarity measure based on gene ontology.} \emph{BMC Bioinformatics} 2016,
#' \strong{17}(1):296)
#' 
#' @importFrom utils setTxtProgressBar txtProgressBar
#' @export
.topo_ic_sim_geneset <- compiler::cmpfun(function(gene_list1,
                                         gene_list2,
                                         topoargs) {
    old <- options(stringsAsFactors = FALSE, warn = -1)
    on.exit(options(old), add = TRUE)
    
    timestart <- Sys.time()
    print(timestart)
    
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
    # apply(unique_pairs, 1, function(pair) {
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
    print(Sys.time()-timestart)
    return(list(GeneSim=score_matrix, 
                GeneList1 = gene_list1, 
                GeneList2 = gene_list2,
                AllGoPairs = topoargs$all_go_pairs))
})

topo_ic_sim_argcheck_genes <- function(ontology, organism, genes1, genes2) {
  # set ontology and organism
  ontology <- match.arg(ontology, c("MF", "BP", "CC"))
  organism <- match.arg(organism, c("human",
                                    "fly",
                                    "mouse",
                                    "rat",
                                    "yeast",
                                    "zebrafish",
                                    "worm",
                                    "arabidopsis",
                                    "ecolik12",
                                    "bovine",
                                    "canine",
                                    "anopheles",
                                    "ecsakai",
                                    "chicken",
                                    "chimp",
                                    "malaria",
                                    "rhesus",
                                    "pig",
                                    "xenopus"))
  if (length(genes1) < 1 || length(genes2) < 1) {
    stop("Not enough genes specified!")
  }
}

topo_ic_sim_genes <- function(ontology,
                              organism,
                              genes1, 
                              genes2, 
                              drop = NULL,
                              progress_bar = T,
                              garbage_collection = F,
                              all_go_pairs = NULL) {
  old <- options(stringsAsFactors = FALSE, warn = -1)
  on.exit(options(old), add = TRUE)
  topo_ic_sim_argcheck_genes(ontology, organism, genes1, genes2)
  # if everything is ok, start preparing...
  topoargs <<- .prepare_variables_topoicsim(organism, ontology, genes1, genes2, 
                                           drop, progress_bar, garbage_collection, all_go_pairs)
  if (length(genes1) == 1 && length(genes2) == 1) {
    # single gene topo
    return(.topo_ic_sim_g1g2(genes1, genes2, topoargs))
  } else {
    # multi gene topo
    return(.topo_ic_sim_geneset(genes1, genes2, topoargs))
  }
}
# 
# topo_ic_sim_term <- function(ontology,
#                              organism,
#                              go1, go2) {
#   topoargs <<- .prepare_variables_topoicsim(organism, ontology, genes1, genes2, 
#                                             drop, progress_bar, garbage_collection, all_go_pairs, )
# }