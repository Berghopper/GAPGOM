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
topo_ic_sim_titj <- compiler::cmpfun(function(go_id1,
                                               go_id2,
                                               ontology,
                                               organism,
                                               weighted_dag,
                                               go_annotation,
                                               root,
                                               IC,
                                               all_info_go_pairs =
                                                 data.frame(GO_1 = 0,
                                                            GO_2 = 0,
                                                            S_1 = 0)) {
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
      all_info_go_pairs <- rbind(all_info_go_pairs, c(go_id1, go_id2, sim))
      return(list(score = sim, go_pairs = all_info_go_pairs))
  } else {
      all_info_go_pairs <- rbind(all_info_go_pairs, c(go_id1, go_id2, 0))
      return(list(score = 0, go_pairs = all_info_go_pairs))
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
                                            go_data,
                                            ontology = "MF",
                                            organism = "yeast",
                                            drop = NULL,
                                            all_info_go_pairs =
                                              data.frame(GO_1 = 0,
                                                         GO_2 = 0,
                                                         S_1 = 0)) {
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

  xx_parents <- switch(ontology, MF = toTable(GOMFPARENTS),
                       BP = toTable(GOBPPARENTS), CC = toTable(GOCCPARENTS))

  go_annotation <- switch(ontology, MF = GOMFANCESTOR, BP = GOBPANCESTOR,
                          CC = GOCCANCESTOR)

  root <- switch(ontology, MF = "GO:0003674", BP = "GO:0008150",
                 CC = "GO:0005575 ")
  weighted_dag <- ftM2graphNEL(as.matrix(xx_parents[, 1:2]))

  goAnno <- go_data@geneAnno
  goAnno <- goAnno[!goAnno$EVIDENCE %in% drop, ]
  go1 <- as.character(unique(goAnno[goAnno[, 1] == gene1, "GO"]))
  go2 <- as.character(unique(goAnno[goAnno[, 1] == gene2, "GO"]))

  if (sum(!is.na(go1)) == 0 || sum(!is.na(go2)) == 0) {
    return(list(geneSim = NA, GO1 = go1, GO2 = go2))
  }
  m <- length(go1)
  n <- length(go2)
  scores <- matrix(nrow = m, ncol = n)
  rownames(scores) <- go1
  colnames(scores) <- go2
  IC <- go_data@IC
  if (gene1 == gene2)
      return(1)
  for (i in 1:m) {
    for (j in 1:n) {
      if (go1[i] == go2[j]) {
        scores[i, j] <- 1
      } else {
        # get the subset of the GO pair.
        DF <- subset(all_info_go_pairs,
                     (all_info_go_pairs[, 1] == go1[i] &
                        all_info_go_pairs[, 2] == go2[j]))
        # now check if they are present, if 0 rows are returned, then check the 
        # inverse combination
        if (nrow(DF) == 0) {
          DF <- subset(all_info_go_pairs,
                       (all_info_go_pairs[, 1] == go2[j] &
                          all_info_go_pairs[, 2] == go1[i]))
        }
        # if this is not returning 0, set score to whatever it already
        # is in the scores df.
        if (nrow(DF) != 0) {
          scores[i, j] <- as.numeric(DF[1, 3])
        } else {
          # if this is not the case (row is not present), then run topo_ic_sim 
          # between 2 terms.
          topo_ic_sim_term_result <-
            topo_ic_sim_titj(go1[i],
                              go2[j],
                              ontology,
                              organism,
                              weighted_dag,
                              go_annotation,
                              root,
                              IC,
                              all_info_go_pairs = all_info_go_pairs)
          scores[i, j] <- topo_ic_sim_term_result$score
          all_info_go_pairs <- topo_ic_sim_term_result$go_pairs
        }
      }
    }
  }
  if (!sum(!is.na(scores)))
      return(list(geneSim = NA, GO1 = go1, GO2 = go2))
  scores <- sqrt(scores)
  sim <- max(sum(sapply(1:m, function(x) {
      max(scores[x, ], na.rm = TRUE)
  }))/m, sum(sapply(1:n, function(x) {
      max(scores[, x], na.rm = TRUE)
  }))/n)
  sim <- round(sim, digits = 3)
  return(sim)
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
    go_data <- set_go_data(organism = organism, ontology = ontology)
    r1 <- length(gene_list1)
    r2 <- length(gene_list2)
    score_matrix <- matrix(nrow = r1, ncol = r2)
    rownames(score_matrix) <- gene_list1
    colnames(score_matrix) <- gene_list2
    for (i in 1:r1) {
        for (j in 1:r2) {
            score_matrix[i, j] <- topo_ic_sim_g1g2(gene_list1[i],
                                                   gene_list2[j], 
                                                   go_data, 
                                                   ontology,
                                                   organism, drop)
        }
    }
    endtime <- Sys.time() - timestart
    cat("took; ", endtime, "s\n")
    return(score_matrix)
})
