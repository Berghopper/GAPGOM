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
  if (is.na(p1) || is.na(p2))
    return(NA)
  if (p1 == 0 || p2 == 0)
    return(NA)

  ancestor1 <- unlist(go_annotation[[go_id1]], FALSE, FALSE)
  ancestor2 <- unlist(go_annotation[[go_id2]], FALSE, FALSE)
  # always intersect all GOIDs, even if they are the same.
  common_ancestor <- intersect(ancestor1, ancestor2)
  return(common_ancestor)
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
  weighted_subgraph <- .set_edge_weight(weighted_subgraph, topoargs$IC)
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


.filter_last_common_ancestors<-function(go1, go2, ontology, organism, root, GA, information_content){
  COMANC <- .common_ancestors(go1, go2, ontology, organism, GA, information_content)
  L=length(COMANC)
  if (length(COMANC)!=0 & !is.na(COMANC)) return(data.frame(GO1=rep(go1, L), GO2=rep(go2, L), LCA=COMANC))
}

.get_short_path_info <- function(GOID, LCA, WDG, GA, information_content){
  sg <- subGraph(c(get(GOID, GA), GOID), WDG)
  sg <- .set_edge_weight(sg, information_content)
  sg <- igraph.to.graphNEL(sg)
  sp <- sp.between(sg, GOID, LCA)
  ic_sp <- sp[[1]]$length
  length_sp <- length(sp[[1]]$path_detail)
  ICGO <- information_content[GOID]
  if (!is.na(ICGO) & ICGO!=0) ic_sp <- ic_sp+(1/(2*ICGO))
  return(data.frame(ICGo_x=ic_sp, LGo_x=length_sp))
}

.all_go_similarities <-function(unique_pairs_genes, topoargs, drop=NULL) {
  organism <- topoargs$organism
  ontology <- topoargs$ontology
  
  GA <- topoargs$go_annotation
  WDG <- topoargs$weighted_dag
  root <- topoargs$root
  # make wa dataframe of all pairs GOs that needed to compute
  All_pairs_GOs <- data.frame(GO1="go2", GO2="go2", stringsAsFactors=F)
  goAnno <- topoargs$go_annotation
  
  all_gos1 <- c()
  all_gos2 <- c()
  for (i in seq_len(nrow(unique_pairs_genes))) {
    pair <- unique_pairs_genes[i,]
    gene1 <- pair[[1]]
    gene2 <- pair[[2]]
    gos1 <- as.character(topoargs$translation_to_goids[
      topoargs$translation_to_goids$ID==gene1,]$GO)
    gos2 <- as.character(topoargs$translation_to_goids[
      topoargs$translation_to_goids$ID==gene2,]$GO)
    all_gos1 <- c(all_gos1, gos1)
    all_gos2 <- c(all_gos2, gos2)
  }
  All_pairs_GOs <- .unique_combos(all_gos1, all_gos2)
  colnames(All_pairs_GOs)<-c("GO1", "GO2")
  
  # find last common ancestors (LCAs)
  ALL_LCAs <- apply(All_pairs_GOs[,c(1,2)], 1, function(x) .filter_last_common_ancestors(x[1], x[2], ontology, organism, root, GA, topoargs$IC))
  # adding LCAs to go pairs
  if(is.null(ALL_LCAs)) Pairs_GOs_LCAs <- as.data.frame(cbind("go1","go2","lca"))
  else Pairs_GOs_LCAs <- as.data.frame(bind_rows(ALL_LCAs))
  
  # remove "all" and "root" as LCAs
  Pairs_GOs_LCAs<-Pairs_GOs_LCAs[(Pairs_GOs_LCAs[,3]!="all"),]
  Pairs_GOs_LCAs<-Pairs_GOs_LCAs[(Pairs_GOs_LCAs[,3]!=root),]
  colnames(Pairs_GOs_LCAs)<-c("GO", "GO", "LCA")
  # convert "go1 go2->LCA" to "go1->LCA and go2->LCA"
  GOs_LCAs <- rbind(Pairs_GOs_LCAs[,c(1,3)], Pairs_GOs_LCAs[,c(2,3)])
  colnames(Pairs_GOs_LCAs)<-c("GO1", "GO2", "LCA")
  GOs_LCAs <- unique(GOs_LCAs)
  
  # compute scores (IC and short pathe length) for "go1->LCA and go2->LCA"
  Pairs_GOs_LCAs_Scores_ <- apply(GOs_LCAs[,c(1,2)], 1, function(x) {
    InfoGO <- .get_short_path_info(x[1], x[2], WDG, GA, topoargs$IC)
    ic_sp <- as.numeric(InfoGO$ICGo_x) 
    length_sp <- as.numeric(InfoGO$LGo_x)
    return(data.frame(GO1=x[1], LCA=x[2], ICSP=ic_sp, LengthSP=length_sp))
  })
  Pairs_GOs_LCAs_Scores <- do.call("rbind", Pairs_GOs_LCAs_Scores_)
  # LCAs
  LCAs<-unique(GOs_LCAs[,2])
  LCAs<-LCAs[!is.na(LCAs)]
  LCAs<-LCAs[LCAs!="all"]
  LCAs<-LCAs[LCAs!=root]
  # compute scores for longest path from LCAs to root
  LCAs_Scores_ <- sapply(LCAs, function(x) {
    wLP_x_root <- .longest_path(topoargs$weighted_dag, topoargs$go_annotation, x, topoargs$root, topoargs$IC)
    return(cbind(LCA=x, wLP=wLP_x_root))
  })
  LCAs_Scores <- as.data.frame(cbind(LCA=LCAs_Scores_[1,], wLP=LCAs_Scores_[2,]),stringsAsFactors=F)
  rownames(LCAs_Scores)<-c()
  
  # merge scores
  AA<-merge(Pairs_GOs_LCAs, Pairs_GOs_LCAs_Scores, by = NULL,by.x = c("GO1", "LCA"), by.y = c("GO1", "LCA"))
  AA1<-merge(AA, Pairs_GOs_LCAs_Scores, by = NULL,by.x = c("GO2", "LCA"), by.y = c("GO1", "LCA"))
  AA2<-merge(AA1, LCAs_Scores, by = NULL,by.x = "LCA", by.y = "LCA")
  
  # compute final scores
  wSP_ti_tj_x <- (AA2$ICSP.x + AA2$ICSP.y)*(AA2$LengthSP.x + AA2$LengthSP.y - 2)
  D_ti_tj <- wSP_ti_tj_x/as.numeric(AA2$wLP)
  
  AA4<-data.frame(GO1=AA2$GO1, GO2=AA2$GO2, D=D_ti_tj)
  AA5 <- summarise(group_by(AA4, GO1, GO2), a_min=min(D))
  Sim_ti_tj <- round(1-(atan(AA5$a_min)/(pi/2)), 3)
  AA5 <- data.frame(AA5, Sim_ti_tj)
  AA5 <- AA5[,-3]
  # AA5 = final dataframe as "GO1 GO2 Sim_ti_tj"
  # update all go pairs
  #assign("AA5", AA5, .GlobalEnv)
  for (i in seq_len(nrow(AA5))) {
    row <- AA5[i,]
    go1 <- row[[1]]
    go2 <- row[[2]]
    topo_score <- row[[3]]
    topoargs$all_go_pairs <- .set_values(go1, go2, topoargs$all_go_pairs, topo_score)
  }
  return(topoargs$all_go_pairs)
}
