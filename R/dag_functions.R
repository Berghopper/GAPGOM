# Directed acyclic graph (DAG) mathematical functions/DAG navigation functions.

# Common Ancestor for two GOIDs
common_ancestor <- function(go_id1, go_id2, ontology, organism, go_annotation, information_content) {
    root_count <- max(information_content[information_content != Inf])
    information_content["all"] = 0
    p1 <- information_content[go_id1]/root_count
    p2 <- information_content[go_id2]/root_count
    if (is.na(p1) || is.na(p2))
        return(NA)
    if (p1 == 0 || p2 == 0)
        return(NA)

    ancestor1 <- unlist(go_annotation[[go_id1]])
    ancestor2 <- unlist(go_annotation[[go_id2]])
    # unlike ppipre, the two GOIDs are always intersected (always the assumption that GOIDs will be different?)  YES, BECAUSE THE SAME GOIDs HAVE PERFECT SIMILARITY 1 BASED ON
    # 'Topoinformation_contentSim_g1g2' FUNCTION.
    common_ancestor <- intersect(ancestor1, ancestor2)
    return(common_ancestor)
}

# Weighting subgraphs to get shortest and longest paths. Weighting edges is done by assigning average inverse half information_content from two corresponding nodes.

# sub_graph = subgraph M = Edge matrix wm = weighted edge matrix

set_edge_weight <- function(sub_graph, information_content) {
    M <- graph::edgeMatrix(sub_graph)
    wm <- graph::eWV(sub_graph, graph::edgeMatrix(sub_graph), useNNames = TRUE)
    Weights <- c()
    for (i in 1:length(names(wm))) {
        go <- strsplit(names(wm)[i], "->")
        ic1 <- information_content[go[[1]][1]]
        ic2 <- information_content[go[[1]][2]]
        w1 <- 0
        w2 <- 0
        if (!is.na(ic1) & ic1 != 0)
            w1 <- 1/(2 * ic1)
        if (!is.na(ic2) & ic2 != 0)
            w2 <- 1/(2 * ic2)
        Weights <- c(Weights, (w1 + w2))
    }
    sub_igraph <- igraph::igraph.from.graphNEL(sub_graph)
    sub_igraph <- igraph::set.edge.attribute(sub_igraph, name = "weight", value = Weights)
    return(sub_igraph)
}

# Calculate weighted long path (winformation_content) between root and a disjunctive common ancestor by topological sorting algorithm.
longest_path <- function(g, lca, root, information_content) {

    tsg <- igraph::topological.sort(g)
    # Set root path attributes Root distance
    igraph::V(g)[tsg[1]]$rdist <- 0
    # Path to root
    igraph::V(g)[tsg[1]]$rpath <- tsg[1]
    # Get longest path from root to every node
    L = 0
    for (node in tsg[-1]) {
        if (tsg[node][[1]]$name == root)
            break
        # Get distance from node's predecessors
        w <- igraph::E(g)[to(node)]$weight
        # Get distance from root to node's predecessors
        d <- igraph::V(g)[nei(node, mode = "in")]$rdist
        # Add distances (assuming one-one corr.)
        wd <- w + d
        # Set node's distance from root to max of added distances
        mwd <- max(wd)
        igraph::V(g)[node]$rdist <- mwd
        # Set node's path from root to path of max of added distances
        mwdn <- as.vec1tor(igraph::V(g)[nei(node, mode = "in")])[match(mwd, wd)]
        V(g)[node]$rpath <- list(c(unlist(igraph::V(g)[mwdn]$rpath), node))
        L <- length(V(g)[node]$rpath[[1]]) - 1
    }
    # Longest path length is the largest distance from root
    lpl <- max(igraph::V(g)$rdist, na.rm = TRUE)
    information_content_lca <- information_content[lca][[1]]
    if (!is.na(information_content_lca) & information_content_lca != 0)
        lpl <- lpl + (1/(2 * information_content_lca))
    return(lpl * L)
}
