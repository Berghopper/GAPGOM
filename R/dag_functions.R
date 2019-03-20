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
#' @return output is different on a case-to-case basis
#'
#' @name dag_funcs
#' @keywords internal
NULL

#' Common Ancestor for two GOIDs
#' @rdname dag_funcs
.common_ancestors <- function(go_id1, go_id2, ontology, organism, go_annotation,
  information_content) {
  root_count <- max(information_content[information_content != Inf])
  information_content["all"] = 0
  p1 <- information_content[go_id1]/root_count
  p2 <- information_content[go_id2]/root_count
  if (is.na(p1) | is.na(p2))
    return(NA)
  if (p1 == 0 | p2 == 0)
    return(NA)

  ancestor1 <- unlist(go_annotation[[go_id1]], FALSE, FALSE)
  ancestor2 <- unlist(go_annotation[[go_id2]], FALSE, FALSE)
  # always intersect all GOIDs, even if they are the same.
  common_ancestors <- intersect(ancestor1, ancestor2)
  return(common_ancestors)
}

#' Weighting subgraphs to get shortest and longest paths. Weighting edges is 
#' done by assigning average inverse half information_content from two 
#' corresponding nodes.
#' @importFrom graph edgeMatrix eWV
#' @importFrom igraph igraph.from.graphNEL set.edge.attribute
#' @rdname dag_funcs
.set_edge_weight <- function(sub_graph, information_content) {
  edge_matrix <- edgeMatrix(sub_graph)
  weighted_edge_matrix <- eWV(sub_graph, edge_matrix, useNNames = TRUE)
  weights <- unlist(
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
    return(w1 + w2)
  }), FALSE, FALSE)
  sub_igraph <- igraph.from.graphNEL(sub_graph)
  sub_igraph <- set.edge.attribute(sub_igraph, 
                                   name = "weight", 
                                   value = weights)
    return(sub_igraph)
}

#' Calculate weighted long path (winformation_content) between root and a 
#' disjunctive common ancestor by topological sorting algorithm.
#' @import igraph
#' @rdname dag_funcs
.longest_path <- function(weighted_dag, go_annotation, last_common_ancestor, 
  root, information_content) {
  # set internal igraph functions to null to appease R CMD check
  to <- NULL
  nei <- NULL
  
  # make weighted subgraph (Subgraph from a disjunctive common ancestor to root)
  weighted_subgraph <- subGraph(c(get(last_common_ancestor, go_annotation), 
    last_common_ancestor), weighted_dag)
  weighted_subgraph <- .set_edge_weight(weighted_subgraph, information_content)
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
}

#' Gets last common ancestors in a special data.frame format.
#' @rdname dag_funcs
.get_last_common_ancestors<-function(go1, go2, ontology, organism, root, 
  go_annotation, information_content) {
  common_ancestors <- .common_ancestors(go1, go2, ontology, organism,
    go_annotation, information_content)
  len_comanc <- length(common_ancestors)
  if (len_comanc != 0 & all(!is.na(common_ancestors))) {
    return(data.table(GO1=rep(go1, len_comanc), GO2=rep(go2, len_comanc), 
      LCA=common_ancestors)) # only return df if there's results. otherwise NA.
  }
}

#' Gets short path info between a go_id and its last common ancestor.
#' @importFrom igraph igraph.from.graphNEL
#' @importFrom graph ftM2graphNEL subGraph
#' @importFrom AnnotationDbi get
#' @importFrom RBGL sp.between
#' @rdname dag_funcs
.get_short_path_info <- function(go_id, last_common_ancestor, weighted_graph, 
  go_annotation, information_content) {
  sg <- subGraph(c(get(go_id, go_annotation), go_id), weighted_graph)
  sg <- .set_edge_weight(sg, information_content)
  sg <- igraph.to.graphNEL(sg)
  sp <- sp.between(sg, go_id, last_common_ancestor)
  ic_sp <- sp[[1]]$length
  length_sp <- length(sp[[1]]$path_detail)
  ic_go <- information_content[go_id]
  if (!is.na(ic_go) & ic_go != 0) {
    # return only if ic is present/not 0.
    ic_sp <- ic_sp+(1/(2*ic_go))
    return(data.frame(ICGo_x=ic_sp, LGo_x=length_sp))
  } 
}

#' Evaluaties all_go_pairs_df and calculates topoicsim scores for only
#' necessary go's. E.g. go pairs that will be 0 will be skipped.
#' @importFrom utils txtProgressBar getTxtProgressBar setTxtProgressBar
#' @importFrom dplyr summarise group_by
#' @importFrom data.table as.data.table rbindlist
#' @rdname dag_funcs
.all_go_similarities <-function(all_go_pairs_df, topoargs, drop=NULL, 
  verbose=FALSE) {
  # set internal df refrences to null to appease R CMD check
  GO1 <- NULL
  GO2 <- NULL
  Distance <- NULL
  
  if (verbose) {message("Started calculating all go's.")}
  
  # Filter out go's present in pre-calculation or ones already present in
  # all_go_pairs (happens if user has specified it).
  if (verbose) {message("Filtering out precalculated values...")}
  removal_rows <- c()
  for (i in seq_len(nrow(all_go_pairs_df))) {
    pair <- all_go_pairs_df[i,]
    go1 <- pair[[1]]
    go2 <- pair[[2]]
    if (!is.na(topoargs$all_go_pairs[go1, go2])) {
      removal_rows <- c(removal_rows, i)
    } else if (!is.na(topoargs$all_go_pairs[go2, go1])) {
      removal_rows <- c(removal_rows, i)
    } else if (topoargs$use_precalculation) {
      if (go1 %in% colnames(topoargs$selected_freq_go_pairs) &
        go2 %in% colnames(topoargs$selected_freq_go_pairs)) {
        # set precalculated value.
        topoargs$all_go_pairs <- .set_values(go1, go2, topoargs$all_go_pairs,
          topoargs$selected_freq_go_pairs[go1, go2])
        # set row for removal.
        removal_rows <- c(removal_rows, i)
      }
    }
  }
  # filter out gos (if necessary)
  if (!is.null(removal_rows)) {
    all_go_pairs_df <- as.data.frame(all_go_pairs_df)[-removal_rows,]
  }
  # set progress bar
  if (topoargs$progress_bar) {
    message("Resolving all common ancestors...")
    pb <- txtProgressBar(min = 0, max = nrow(all_go_pairs_df), style = 3)
  }
  # find last common ancestors (lcas)
  all_lcas <- apply(all_go_pairs_df[, c(1, 2)], 1, function(go_pair) {
    if (topoargs$progress_bar) {
      setTxtProgressBar(pb, getTxtProgressBar(pb)+1)
    }
    return(.get_last_common_ancestors(go_pair[1], go_pair[2], topoargs$ontology, 
       topoargs$organism, topoargs$root, topoargs$go_annotation, topoargs$IC))
    })
  # adding lcas to go pairs
  if (is.null(all_lcas)) {
    go_lca_pairs <- as.data.frame(cbind("GO1", "GO2", "LCA"))  
  }
  else {
    go_lca_pairs <- as.data.frame(rbindlist(all_lcas)) 
  }
  
  # remove "all" and "root" as lca's
  go_lca_pairs <- go_lca_pairs[(go_lca_pairs[,3] != "all"),]
  go_lca_pairs <- go_lca_pairs[(go_lca_pairs[,3] != topoargs$root),]
  colnames(go_lca_pairs)<-c("GO", "GO", "LCA")
  # convert "go1 go2 -> LCA" to "go1->LCA and go2->LCA" so we get all unique
  # go - lca pairs
  gos_lcas <- rbind(go_lca_pairs[,c(1, 3)], go_lca_pairs[,c(2, 3)])
  colnames(go_lca_pairs)<-c("GO1", "GO2", "LCA")
  gos_lcas <- unique(gos_lcas)
  
  # prepare lca's
  lcas <- unique(gos_lcas[,2])
  lcas <- lcas[!is.na(lcas)]
  lcas <- lcas[lcas!="all"]
  lcas <- lcas[lcas!=topoargs$root]
  
  if (topoargs$progress_bar) {
    message("\nDone!")
    message("Calculating short paths...")
    pb <- txtProgressBar(min = 0, max = nrow(gos_lcas), style = 3)
  }
  # compute scores (IC and short path length) for "go1->LCA and go2->LCA"
  go_lca_pair_scores_list <- apply(gos_lcas[,c(1, 2)], 1, function(go_terms) {
    InfoGO <- .get_short_path_info(go_terms[1], go_terms[2], 
      topoargs$weighted_dag, topoargs$go_annotation, topoargs$IC)
    ic_sp <- as.numeric(InfoGO$ICGo_x) 
    length_sp <- as.numeric(InfoGO$LGo_x)
    if (topoargs$progress_bar) {
      setTxtProgressBar(pb, getTxtProgressBar(pb)+1)
    }
    return(data.table(GO1=go_terms[1], LCA=go_terms[2], ICSP=ic_sp, 
      LengthSP=length_sp))
  })
  # bind results
  go_lca_pair_scores <- rbindlist(go_lca_pair_scores_list)
  
  if (topoargs$progress_bar) {
    message("\nDone!")
    message("Calculating long paths...")
    pb <- txtProgressBar(min = 0, max = length(lcas), style = 3)
  }
  # compute scores for longest path from lcas to root
  lca_scores_list <- sapply(lcas, function(lca) {
    wLP_x_root <- .longest_path(topoargs$weighted_dag, topoargs$go_annotation,
      lca, topoargs$root, topoargs$IC)
    if (topoargs$progress_bar) {
      setTxtProgressBar(pb, getTxtProgressBar(pb)+1)
    }
    return(cbind(LCA=lca, wLP=wLP_x_root))
  })
  lca_scores <- as.data.frame(cbind(LCA=lca_scores_list[1,], 
                                    wLP=lca_scores_list[2,]))
  rownames(lca_scores) <- c()
  
  if (topoargs$progress_bar) {message("\nMerging...")}
  # merge scores
  merged_scores <- merge(go_lca_pairs, go_lca_pair_scores, by = NULL, 
    by.x = c("GO1", "LCA"), by.y = c("GO1", "LCA"))
  merged_scores <- merge(merged_scores, go_lca_pair_scores, by = NULL, 
    by.x = c("GO2", "LCA"), by.y = c("GO1", "LCA"))
  merged_scores <- merge(merged_scores, lca_scores, by = NULL, by.x = "LCA", 
    by.y = "LCA")
  
  # compute final scores
  wSP_ti_tj_x <- (merged_scores$ICSP.x + merged_scores$ICSP.y)*
    (merged_scores$LengthSP.x + merged_scores$LengthSP.y - 2)
  D_ti_tj <- wSP_ti_tj_x/as.numeric(merged_scores$wLP)
  
  go_distance_df <- data.frame(GO1=merged_scores$GO1, GO2=merged_scores$GO2,
    Distance=D_ti_tj)
  go_similarity_df <- summarise(group_by(go_distance_df, GO1, GO2), 
    a_min=min(Distance))
  sim_ti_tj <- round(1-(atan(go_similarity_df$a_min)/(pi/2)), 3)
  go_similarity_df <- data.frame(go_similarity_df, sim_ti_tj)
  go_similarity_df <- go_similarity_df[,-3]
  
  # AA5 = final dataframe as "GO1 GO2 Sim_ti_tj"
  # update all go pairs
  for (i in seq_len(nrow(go_similarity_df))) {
    row <- go_similarity_df[i,]
    go1 <- row[[1]]
    go2 <- row[[2]]
    topo_score <- row[[3]]
    topoargs$all_go_pairs <- .set_values(go1, go2, topoargs$all_go_pairs, 
      topo_score)
  }
  if (verbose) {message("Done!")}
  return(topoargs$all_go_pairs)
}
