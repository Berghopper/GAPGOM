#' topo_ic_sim algorithm to calculate similarity between two GO terms.
#' @importFrom GO.db GOMFCHILDREN GOBPCHILDREN GOCCCHILDREN GOMFPARENTS
#' GOBPPARENTS GOCCPARENTS GOMFANCESTOR GOBPANCESTOR GOCCANCESTOR
#' @importFrom igraph igraph.from.graphNEL
#' @importFrom graph ftM2graphNEL subGraph
#' @importFrom AnnotationDbi get
#' @importFrom RBGL sp.between
#' @export
topo_ic_sim_ti_tj <- function(go_id1, go_id2, ont, organism, weighted_dag, go_annotation, root, IC, all_info_go_pairs = data.frame(GO_1 = 0, GO_2 = 0, S_1 = 0)) {
    # =data.frame(GO_1=0, GO_2=0, S_1=0)
    old <- options(stringsAsFactors = FALSE, warn = -1)
    on.exit(options(old), add = TRUE)
    D_ti_tj_x <- c()
    common_ancestor <- common_ancestor(go_id1, go_id2, ont, organism, go_annotation, IC)
    if (length(common_ancestor) != 0 && !is.na(common_ancestor)) {
        for (x in common_ancestor) {
            # To identify all disjunctive common ancestors
            immediate_children_x <- switch(ont, MF = GOMFCHILDREN[[x]], BP = GOBPCHILDREN[[x]], CC = GOCCCHILDREN[[x]])
            if (x != "all" & x != root & !is.na(x) & length(intersect(immediate_children_x, common_ancestor)) == 0) {
                # Subgraph from two GO terms go_id1 and go_id2
                sg1 <- subGraph(c(get(go_id1, go_annotation), go_id1), weighted_dag)
                sg2 <- subGraph(c(get(go_id2, go_annotation), go_id2), weighted_dag)
                # Subgraph from a disjunctive common ancestor to root
                sglca <- subGraph(c(get(x, go_annotation), x), weighted_dag)
                sglca <- set_edge_weight(sglca, IC)
                wLP_x_root <- longest_path(sglca, x, root, IC)
                sg1 <- set_edge_weight(sg1, IC)
                sg1 <- igraph.to.graphNEL(sg1)
                sp1 <- sp.between(sg1, go_id1, x)
                ic_sp1 <- sp1[[1]]$length
                length_sp1 <- length(sp1[[1]]$path_detail)

                sg2 <- set_edge_weight(sg2, IC)
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
}

#' topo_ic_sim algorithm to calculate similarity between two genes or gene products.
#' @importFrom AnnotationDbi toTable
#' @importFrom graph ftM2graphNEL
#' @export
topo_ic_sim_g1g2 <- function(gene1, gene2, go_data, ont = "MF", organism = "yeast", drop = NULL, all_info_go_pairs = data.frame(GO_1 = 0, GO_2 = 0, S_1 = 0)) {
    old <- options(stringsAsFactors = FALSE, warn = -1)
    on.exit(options(old), add = TRUE)
    # set ontology and organism
    ont <- match.arg(ont, c("MF", "BP", "CC"))
    organism <- match.arg(organism, c("human", "fly", "mouse", "rat", "yeast", "zebrafish", "worm", "arabidopsis", "ecolik12", "bovine", "canine", "anopheles", "ecsakai", "chicken",
        "chimp", "malaria", "rhesus", "pig", "xenopus"))
    # Initial definitions and sets based on organism and ontology type

    xx_parents <- switch(ont, MF = toTable(GOMFPARENTS), BP = toTable(GOBPPARENTS), CC = toTable(GOCCPARENTS))

    go_annotation <- switch(ont, MF = GOMFANCESTOR, BP = GOBPANCESTOR, CC = GOCCANCESTOR)

    root <- switch(ont, MF = "GO:0003674", BP = "GO:0008150", CC = "GO:0005575 ")
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
            if (go1[i] == go2[j])
                scores[i, j] <- 1 else {
                # get the subset of the GO pair.
                DF <- subset(all_info_go_pairs, (all_info_go_pairs[, 1] == go1[i] & all_info_go_pairs[, 2] == go2[j]))
                # now check if they are present, if 0 rows are returned, then check the inverse combination
                if (nrow(DF) == 0) {
                  DF <- subset(all_info_go_pairs, (all_info_go_pairs[, 1] == go2[j] & all_info_go_pairs[, 2] == go1[i]))
                }
                # if this is not returning 0, set score to whatever it already is in the scores df.
                if (nrow(DF) != 0) {
                  scores[i, j] <- as.numeric(DF[1, 3])
                } else {
                  # if this is not the case (row is not present), then run topo_ic_sim between 2 terms.
                  topo_ic_sim_term_result <- topo_ic_sim_ti_tj(go1[i], go2[j], ont, organism, weighted_dag, go_annotation, root, IC, all_info_go_pairs = all_info_go_pairs)
                  scores[i, j] <- topo_ic_sim_term_result$score
                  # print(scores)
                  all_info_go_pairs <- topo_ic_sim_term_result$go_pairs
                }
            }
        }
    }
    if (!sum(!is.na(scores)))
        return(list(geneSim = NA, GO1 = go1, GO2 = go2))
    scores <- sqrt(scores)  # Some scores seem unnecesary to calculate? Like, when they have the same GOID the similarity will always be 1.
    ## YES, I ADDED A NEW LINE FOR THIS CASE ABOVE //if(gene1==gene2) return(1)// View(scores)
    sim <- max(sum(sapply(1:m, function(x) {
        max(scores[x, ], na.rm = TRUE)
    }))/m, sum(sapply(1:n, function(x) {
        max(scores[, x], na.rm = TRUE)
    }))/n)
    sim <- round(sim, digits = 3)
    return(sim)
    # return (list(geneSim=sim, GO1=go1, GO2=go2))
}


#' Calculate topo_ic_sim between two gene lists
#' @export
topo_ic_sim <- function(gene_list1, gene_list2, ont = "MF", organism = "human", drop = NULL) {
    old <- options(stringsAsFactors = FALSE, warn = -1)
    on.exit(options(old), add = TRUE)
    timestart <- Sys.time()
    go_data <- set_go_data(organism = organism, ontology = ont)
    r1 <- length(gene_list1)
    r2 <- length(gene_list2)
    score_matrix <- matrix(nrow = r1, ncol = r2)
    rownames(score_matrix) <- gene_list1
    colnames(score_matrix) <- gene_list2
    for (i in 1:r1) {
        for (j in 1:r2) {
            score_matrix[i, j] <- topo_ic_sim_g1g2(gene_list1[i], gene_list2[j], go_data, ont, organism, drop)
        }
    }
    endtime <- Sys.time() - timestart
    cat("took; ", endtime, "s\n")
    return(score_matrix)
}
