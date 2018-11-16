#' GAPGOM internal - Directed acyclic graph (DAG) mathematical functions/DAG 
#' navigation functions.
#' 
#' These functions are internal functions and should not be called by the user.
#' 
#' Mathematical calculations and navigations functions for directed acyclic 
#' graphs (DAG). These functions do a multitude of things; calculate longest
#' paths between two nodes, calculate edge weight and retrieving a common 
#' ancestor in the GO DAG given two GO ids.
#' 
#' @section Notes:
#' Internal functions all used in topo_ic_sim_titj().
#'
#' @name dag_funcs
#' @keywords internal
NULL

#' Common Ancestor for two GOIDs
#' @rdname dag_funcs
.common_ancestors <- compiler::cmpfun(function(go_id1, go_id2, ontology, 
                                             organism, go_annotation, 
                                             information_content) {
  root_count <- max(information_content[information_content != Inf])
  information_content["all"] = 0
  p1 <- information_content[go_id1]/root_count
  p2 <- information_content[go_id2]/root_count
  if (is.na(p1) || is.na(p2))
    return(NA)
  if (p1 == 0 || p2 == 0)
    return(NA)

  ancestor1 <- unlist(go_annotation[[go_id1]], FALSE, FALSE)
  ancestor2 <- unlist(go_annotation[[go_id2]], FALSE, FALSE)
  # always intersect all GOIDs, even if they are the same.
  common_ancestor <- intersect(ancestor1, ancestor2)
  return(common_ancestor)
})

#' Weighting subgraphs to get shortest and longest paths. Weighting edges is 
#' done by assigning average inverse half information_content from two 
#' corresponding nodes.
#' @importFrom graph edgeMatrix eWV
#' @importFrom igraph igraph.from.graphNEL set.edge.attribute
#' @rdname dag_funcs
.set_edge_weight <- compiler::cmpfun(function(sub_graph, information_content) {
  edge_matrix <- edgeMatrix(sub_graph)
  weighted_edge_matrix <- eWV(sub_graph, edge_matrix, useNNames = TRUE)
  Weights <- c()
  lapply(seq_along(names(weighted_edge_matrix)), function(i) {
    go <- strsplit(names(weighted_edge_matrix)[i], "->")
    ic1 <- information_content[go[[1]][1]]
    ic2 <- information_content[go[[1]][2]]
    w1 <- 0
    w2 <- 0
    if (!is.na(ic1) & ic1 != 0)
      w1 <- 1/(2 * ic1)
    if (!is.na(ic2) & ic2 != 0)
      w2 <- 1/(2 * ic2)
    Weights <<- c(Weights, (w1 + w2))
  })
  sub_igraph <- igraph.from.graphNEL(sub_graph)
  sub_igraph <- set.edge.attribute(sub_igraph, 
                                   name = "weight", 
                                   value = Weights)
    return(sub_igraph)
})

#' Calculate weighted long path (winformation_content) between root and a 
#' disjunctive common ancestor by topological sorting algorithm.
#' @import igraph
#' @rdname dag_funcs
.longest_path <- compiler::cmpfun(function(weighted_subgraph, 
                                          last_common_ancestor, 
                                          root, 
                                          information_content) {
  tsg <- topological.sort(weighted_subgraph)
  # Set root path attributes Root distance
  V(weighted_subgraph)[tsg[1]]$rdist <- 0
  # Path to root
  V(weighted_subgraph)[tsg[1]]$rpath <- tsg[1]
  # Get longest path from root to every node
  L = 0
  for (node in tsg[-1]) {
    if (tsg[node][[1]]$name == root) {
        break
      }
    # Get distance from node's predecessors
    w <- E(weighted_subgraph)[to(node)]$weight
    # Get distance from root to node's predecessors
    d <- V(weighted_subgraph)[nei(node, mode = "in")]$rdist
    # Add distances (assuming one-one corr.)
    wd <- w + d
    # Set node's distance from root to max of added distances
    mwd <- max(wd)
    V(weighted_subgraph)[node]$rdist <- mwd
    # Set node's path from root to path of max of added distances
    mwdn <- as.vector(V(weighted_subgraph)[nei(node, mode = "in")])[
      match(mwd, wd)]
    V(weighted_subgraph)[node]$rpath <- list(c(unlist(V(weighted_subgraph)[
      mwdn]$rpath), node))
    L <- length(V(weighted_subgraph)[node]$rpath[[1]]) - 1
  }
  # Longest path length is the largest distance from root
  lpl <- max(V(weighted_subgraph)$rdist, na.rm = TRUE)
  information_content_lca <- information_content[last_common_ancestor][[1]]
  if (!is.na(information_content_lca) & information_content_lca != 0)
    lpl <- lpl + (1/(2 * information_content_lca))
  return(lpl * L)
})
