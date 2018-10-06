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
#' @param IC Information content of the two terms
#' @param all_info_go_pairs dataframe of GO Term pairs with a column
#' representing similarity between the two. it is recommended to always keep
#' this on default unless you know what you are doing.
#' @return list of score and all_info_go_pairs, the latter of which will be 
#' reused.
#' @examples
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
                                               ontology,
                                               organism,
                                               weighted_dag,
                                               go_annotation,
                                               root,
                                               IC) {
  old <- options(stringsAsFactors = FALSE, warn = -1)
  on.exit(options(old), add = TRUE)
  D_ti_tj_x <- c()
  common_ancestor <- .common_ancestor(go_id1, go_id2, ontology, organism,
                                     go_annotation, IC)
  if (length(common_ancestor) != 0 && !is.na(common_ancestor)) {
    for (x in common_ancestor) {
      # To identify all disjunctive common ancestors
      immediate_children_x <- switch(ontology,
                                     MF = GOMFCHILDREN[[x]],
                                     BP = GOBPCHILDREN[[x]],
                                     CC = GOCCCHILDREN[[x]])
      
      if (x != "all" & x != root & !is.na(x) &
        length(intersect(immediate_children_x, common_ancestor)) == 0) {
        # Subgraph from two GO terms go_id1 and go_id2
        sg1 <- subGraph(c(get(go_id1, go_annotation), go_id1),
                        weighted_dag)
        sg2 <- subGraph(c(get(go_id2, go_annotation), go_id2),
                        weighted_dag)
        # Subgraph from a disjunctive common ancestor to root
        sglca <- subGraph(c(get(x, go_annotation), x), weighted_dag)
        sglca <- .set_edge_weight(sglca, IC)
        wLP_x_root <- .longest_path(sglca, x, root, IC)
        sg1 <- .set_edge_weight(sg1, IC)
        sg1 <- igraph.to.graphNEL(sg1)
        sp1 <- sp.between(sg1, go_id1, x)
        ic_sp1 <- sp1[[1]]$length
        length_sp1 <- length(sp1[[1]]$path_detail)
        sg2 <- .set_edge_weight(sg2, IC)
        sg2 <- igraph.to.graphNEL(sg2)
        sp2 <- sp.between(sg2, go_id2, x)
        ic_sp2 <- sp2[[1]]$length
        length_sp2 <- length(sp2[[1]]$path_detail)
        IC_GOID_1 <- IC[go_id1][[1]]
        IC_GOID_2 <- IC[go_id2][[1]]
        if (!is.na(IC_GOID_1) & IC_GOID_1 != 0)
          ic_sp1 <- ic_sp1 + (1/(2 * IC_GOID_1))
        if (!is.na(IC_GOID_2) & IC_GOID_2 != 0)
          ic_sp2 <- ic_sp2 + (1/(2 * IC_GOID_2))

        wSP_ti_tj_x <- (ic_sp1 + ic_sp2) * (length_sp1 + length_sp2 - 2)
        D_ti_tj_x <- c(D_ti_tj_x, wSP_ti_tj_x/wLP_x_root)
      }
    }
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
#' @param all_info_go_pairs dataframe of GO Term pairs with a column
#' representing similarity between the two. it is recommended to always keep
#' this on default unless you know what you are doing.
#'
#' @return similarity between genes taken from the mean of all term 
#' similarities.
#' @examples
#' GAPGOM::topo_ic_sim_g1g2("218", "501", ont="MF", organism="human", drop=NULL)
#'
#' @references [1] Ehsani R, Drablos F: \strong{TopoICSim: a new semantic
#' similarity measure based on gene ontology.} \emph{BMC Bioinformatics} 2016,
#' \strong{17}(1):296)
#' @importFrom AnnotationDbi toTable
#' @importFrom graph ftM2graphNEL
#' @export
topo_ic_sim_g1g2 <- compiler::cmpfun(function(gene1,
                                            gene2,
                                            ontology = "MF",
                                            organism = "yeast",
                                            go_data = NULL,
                                            drop = NULL,
                                            translation_to_goids = NULL,
                                            all_info_go_pairs = NULL
                                            ) {
  old <- options(stringsAsFactors = FALSE, warn = -1)
  on.exit(options(old), add = TRUE)
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
  # Initial definitions and sets based on organism and ontology type
  if (is.null(go_data)) {
    go_data <- .set_go_data(organism = organism, ontology = ontology) 
  }
  
  xx_parents <- switch(ontology, MF = toTable(GOMFPARENTS),
                       BP = toTable(GOBPPARENTS), CC = toTable(GOCCPARENTS))

  go_annotation <- switch(ontology, MF = GOMFANCESTOR, BP = GOBPANCESTOR,
                          CC = GOCCANCESTOR)

  root <- switch(ontology, MF = "GO:0003674", BP = "GO:0008150",
                 CC = "GO:0005575 ")
  weighted_dag <- ftM2graphNEL(as.matrix(xx_parents[, 1:2]))
  
  if (is.null(translation_to_goids)) {
    # lookup goids of all gene ids
    translation_to_goids <- .go_ids_lookup(unique(c(gene_list1, gene_list2)), 
                                           go_data, 
                                           drop = drop)
  }
  if (is.null(all_info_go_pairs)) {
    go_unique_list <- unique(translation_to_goids$GO)
    all_info_go_pairs <- .prepare_score_matrix_topoicsim(go_unique_list, 
                                                         go_unique_list)
  }
  gos1 <<- as.character(translation_to_goids[translation_to_goids$ID==gene1,]$GO)
  gos2 <<- as.character(translation_to_goids[translation_to_goids$ID==gene2,]$GO)

  if (sum(!is.na(gos1)) == 0 || sum(!is.na(gos2)) == 0) {
    return(list(geneSim = NA, GO1 = gos1, GO2 = gos2, all_info_go_pairs = all_info_go_pairs))
  }
  scores <- .prepare_score_matrix_topoicsim(gos1, gos2)
  IC <- go_data@IC
  
  unique_pairs <<- .unique_combos(gos1, gos2)
  
  #View(all_info_go_pairs)
  #View(unique_pairs)
  apply(unique_pairs, 1, function(pair) {
    go1 <- pair[1]
    go2 <- pair[2]
    # if this is not the case (row is not present), then run topo_ic_sim 
    # between 2 terms.
    if (is.na(all_info_go_pairs[go1, go2])) {
      score <-
        .topo_ic_sim_titj(go1,
                          go2,
                          ontology,
                          organism,
                          weighted_dag,
                          go_annotation,
                          root,
                          IC)
      #View(scores)
      scores <<- .set_values(go1, go2, scores, score)
      all_info_go_pairs[go1, go2] <<- score
      all_info_go_pairs[go2, go1] <<- score
    } else {
      # set already existing value.
      scores <<- .set_values(go1, go2, scores, all_info_go_pairs[go1, go2])
    }
    
    #print(dim(all_info_go_pairs))
  })
  #print(scores)
  #print(dim(scores))
  #print(unique_pairs)
  if (!sum(!is.na(scores))) {
    #print("uhoh...")
    return(list(geneSim = NA, GO1 = gos1, GO2 = gos2, all_info_go_pairs = all_info_go_pairs))
  }
  scores <- sqrt(scores)
  m <- length(gos1)
  n <- length(gos2)
  sim <- max(sum(sapply(1:m, function(x) {
      max(scores[x, ], na.rm = TRUE)
  }))/m, sum(sapply(1:n, function(x) {
      max(scores[, x], na.rm = TRUE)
  }))/n)
  sim <- round(sim, digits = 3)
  #print(sim)
  return(list(geneSim = sim, GO1 = gos1, GO2 = gos2, all_info_go_pairs = all_info_go_pairs))
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
#' @return returns dataframe containing similarities between corresponding
#' Gene pairs.
#' @examples
#' list1 <- c("126133","221","218","216","8854","220","219","160428","224",
#' "222","8659","501","64577","223","217","4329","10840","7915")
#' GAPGOM::topo_ic_sim(list1, list1, ont="MF", organism="human", drop=NULL)
#' 
#' @references [1] Ehsani R, Drablos F: \strong{TopoICSim: a new semantic
#' similarity measure based on gene ontology.} \emph{BMC Bioinformatics} 2016,
#' \strong{17}(1):296)
#' @export
topo_ic_sim <- compiler::cmpfun(function(gene_list1,
                                         gene_list2,
                                         ontology = "MF",
                                         organism = "human",
                                         drop = NULL) {
    old <- options(stringsAsFactors = FALSE, warn = -1)
    on.exit(options(old), add = TRUE)
    
    timestart <- Sys.time()
    print(timestart)
    
    go_data <- .set_go_data(organism = organism, ontology = ontology)
    # lookup goids of all gene ids
    translation_to_goids <- .go_ids_lookup(unique(c(gene_list1, gene_list2)), 
                                           go_data, 
                                           drop = drop)
    
    # set up score matrix in advance
    score_matrix <- .prepare_score_matrix_topoicsim(gene_list1, gene_list2)
    
    go_unique_list <- unique(translation_to_goids$GO)
    all_info_go_pairs <- .prepare_score_matrix_topoicsim(go_unique_list, 
                                                         go_unique_list)
    
    
    # only loop through necesary vectors (unique pairs)
    unique_pairs <- .unique_combos(gene_list1, gene_list2)
    apply(unique_pairs, 1, function(pair) {
      gene1 <- pair[1]
      gene2 <- pair[2]
      genepair_result <- topo_ic_sim_g1g2(gene1,
                                          gene2, 
                                          ontology,
                                          organism, 
                                          go_data = go_data, 
                                          drop = drop, 
                                          translation_to_goids = 
                                            translation_to_goids,
                                          all_info_go_pairs = all_info_go_pairs)
      
      score_matrix <<- .set_values(gene1, gene2, score_matrix, 
                                  genepair_result$geneSim)
      all_info_go_pairs <<- genepair_result$all_info_go_pairs
    })
    print(Sys.time()-timestart)
    return(score_matrix)
})

