# Directed acyclic graph (DAG) mathematical functions/DAG navigation functions.

# Common Ancestor for two GOIDs
ComAnc <- function(GOID1,GOID2,ont, organism, GA){
  rootCount <- max(IC[IC != Inf])
  IC["all"] = 0
  p1 <- IC[GOID1]/rootCount
  p2 <- IC[GOID2]/rootCount
  if(is.na(p1) || is.na(p2)) return (NA)
  if (p1 == 0 || p2 == 0) return (NA)

  ancestor1 <- unlist(GA[[GOID1]])
  ancestor2 <- unlist(GA[[GOID2]])
  # unlike ppipre, the two GOIDs are always intersected (always the assumption that GOIDs will be different?)
  ## YES, BECAUSE THE SAME GOIDs HAVE PERFECT SIMILARITY 1 BASED ON "TopoICSim_g1g2" FUNCTION.
  commonAncestor <- intersect(ancestor1, ancestor2)
  return(commonAncestor)
}

# Weighting subgraphs to get shortest and longest paths. Weighting edges is done by
# assigning average inverse half IC from two corresponding nodes.

# SG = subgraph
# M = Edge matrix
# wm = weighted edge matrix

WSG<-function(SG){
  M <- edgeMatrix(SG)
  wm <- eWV(SG, edgeMatrix(SG), useNNames=TRUE)
  Weights <- c()
  for (i in 1:length(names(wm))) {
    Go_<- strsplit(names(wm)[i], "->")
    ic1 <- IC[Go_[[1]][1]]
    ic2 <- IC[Go_[[1]][2]]
    w1 <- 0
    w2 <- 0
    if (!is.na(ic1) & ic1!=0)  w1<-1/(2*ic1)
    if (!is.na(ic2) & ic2!=0)  w2<-1/(2*ic2)
    Weights <- c(Weights, (w1+w2))
  }
  SG1 <- igraph.from.graphNEL(SG)
  SG1 <- set.edge.attribute(SG1, name="weight", value=Weights)
  return(SG1)
}

# Calculate weighted long path (wIC) between root and a disjunctive common ancestor
# by topological sorting algorithm.
Longest_Path<-function(g, lca, root){

  tsg <- topological.sort(g)
  # Set root path attributes
  # Root distance
  V(g)[tsg[1]]$rdist <- 0
  # Path to root
  V(g)[tsg[1]]$rpath <- tsg[1]
  # Get longest path from root to every node
  L=0
  for(node in tsg[-1])
  {
    if (tsg[node][[1]]$name==root) break
    # Get distance from node's predecessors
    w <- E(g)[to(node)]$weight
    # Get distance from root to node's predecessors
    d <- V(g)[nei(node,mode="in")]$rdist
    # Add distances (assuming one-one corr.)
    wd <- w+d
    # Set node's distance from root to max of added distances
    mwd <- max(wd)
    V(g)[node]$rdist <- mwd
    # Set node's path from root to path of max of added distances
    mwdn <- as.vector(V(g)[nei(node,mode="in")])[match(mwd,wd)]
    V(g)[node]$rpath <- list(c(unlist(V(g)[mwdn]$rpath), node))
    L<-length(V(g)[node]$rpath[[1]])-1
  }
  # Longest path length is the largest distance from root
  lpl <- max(V(g)$rdist, na.rm=TRUE)
  IC_lca<-IC[lca][[1]]
  if (!is.na(IC_lca) & IC_lca!=0) lpl<-lpl+(1/(2*IC_lca))
  return(lpl * L)
}
