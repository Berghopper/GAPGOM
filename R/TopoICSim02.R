
# R code to measure similarity two genes based on TopoICSim algorithm.
# Ehsani R, Drablos F. TopoICSim: A new semantic similarity measure based on gene ontology. In revision.
# Written by Rezvan Ehsani, rezvan.ehsani@ntnu.no.

###################
#rm(list = ls())

# load libs
library("GO.db")
library("graph")
library("RBGL")
library("igraph")
library("Rgraphviz")
library("ontologySimilarity")
library("ontologyIndex")

options(warn=-1)
options(stringsAsFactors = FALSE)

# Re-used from ppiPre package
# (loads in correct package namespace for species database)
CheckAnnotationPackage <- function(species){
  if (species == "human")
    if(!requireNamespace("org.Hs.eg.db"))
      stop("The package org.Hs.eg.db is needed.")
  if (species == "yeast")
    if(!requireNamespace("org.Sc.sgd.db"))
      stop("The package org.Sc.sgd.db is needed.")
  if (species == "fly")
    if(!requireNamespace("org.Dm.eg.db"))
      stop("The package org.Dm.eg.db is needed.")
  if (species == "mouse")
    if(!requireNamespace("org.Mm.eg.db"))
      stop("The package org.Mm.eg.db is needed.")
  if (species == "rat")
    if(!requireNamespace("org.Rn.eg.db"))
      stop("The package org.Rn.eg.db is needed.")
  if (species == "zebrafish")
    if(!requireNamespace("org.Dr.eg.db"))
      stop("The package org.Dr.eg.db is needed.")
  if (species == "worm")
    if(!requireNamespace("org.Ce.eg.db"))
      stop("The package org.Ce.eg.db is needed.")
  if (species == "arabidopsis")
    if(!requireNamespace("org.At.tair.db"))
      stop("The package org.At.tair.db is needed.")
  if (species == "ecolik12")
    if(!requireNamespace("org.EcK12.eg.db"))
      stop("The package org.EcK12.eg.db is needed.")
  if (species == "bovine")
    if(!requireNamespace("org.Bt.eg.db"))
      stop("The package org.Bt.eg.db is needed.")
  if (species == "canine")
    if(!requireNamespace("org.Cf.eg.db"))
      stop("The package org.Cf.eg.db is needed.")
  if (species == "anopheles")
    if(!requireNamespace("org.Ag.eg.db"))
      stop("The package org.Ag.eg.db is needed.")
  if (species == "ecsakai")
    if(!requireNamespace("org.EcSakai.eg.db"))
      stop("The package org.EcSakai.eg.db is needed.")
  if (species == "chicken")
    if(!requireNamespace("org.Gg.eg.db"))
      stop("The package org.Gg.eg.db is needed.")
  if (species == "chimp")
    if(!requireNamespace("org.Pt.eg.db"))
      stop("The package org.Pt.eg.db is needed.")
  if (species == "malaria")
    if(!requireNamespace("org.Pf.plasmo.db"))
      stop("The package org.Pf.plasmo.db is needed.")
  if (species == "rhesus")
    if(!requireNamespace("org.Mmu.eg.db"))
      stop("The package org.Mmu.eg.db is needed.")
  if (species == "pig")
    if(!requireNamespace("org.Ss.eg.db"))
      stop("The package org.Ss.eg.db is needed.")
  if (species == "xenopus")
    if(!requireNamespace("org.Xl.eg.db"))
      stop("The package org.Xl.eg.db is needed.")
}

# Re-used from ppiPre package
# Get gene ontology measures given a certain ontology of an organism.
GetOntology <-  function(gene, organism, ontology, dropCodes) {
    species <- switch(organism,
                    human = "Hs",
                    fly = "Dm",
                    mouse = "Mm",
                    rat = "Rn",
                    yeast = "Sc",
                    zebrafish = "Dr",
                    worm = "Ce",
                    arabidopsis = "At",
                    ecolik12 = "EcK12",
                    bovine	= "Bt",
                    canine	= "Cf",
                    anopheles	=	"Ag",
                    ecsakai	=	"EcSakai",
                    chicken	=	"Gg",
                    chimp	=	"Pt",
                    malaria	=	"Pf",
                    rhesus	=	"Mmu",
                    pig	= "Ss",
                    xenopus	=	"Xl"
    )

    # grab correct gene ontology map
    gomap <- switch(organism,
                    human = org.Hs.eg.db::org.Hs.egGO,
                    fly = org.Dm.eg.db::org.Dm.egGO,
                    mouse = org.Mm.eg.db::org.Mm.egGO,
                    rat = org.Rn.eg.db::org.Rn.egGO,
                    yeast = org.Sc.sgd.db::org.Sc.sgdGO,
                    zebrafish = org.Dr.eg.db::org.Dr.egGO,
                    worm = org.Ce.eg.db::org.Ce.egGO,
                    arabidopsis = org.At.tair.db::org.At.tairGO,
                    ecoli = org.EcK12.eg.db::org.EcK12.egGO,
                    bovine	= org.Bt.eg.db::org.Bt.egGO,
                    canine	= org.Cf.eg.db::org.Cf.egGO,
                    anopheles	=	org.Ag.eg.db::org.Ag.egGO,
                    ecsakai	=	org.EcSakai.eg.db::org.EcSakai.egGO,
                    chicken	=	org.Gg.eg.db::org.Gg.egGO,
                    chimp	=	org.Pt.eg.db::org.Pt.egGO,
                    malaria	=	org.Pf.plasmo.db::org.Pf.plasmoGO,
                    rhesus	=	org.Mmu.eg.db::org.Mmu.egGO,
                    pig	= org.Ss.eg.db::org.Ss.egGO,
                    xenopus	=	org.Xl.eg.db::org.Xl.egGO
        )
    # check for your specific gene in the database
    allGO <- gomap[[gene]]
    if (is.null(allGO)) {
        return (NA)
        }
    if (sum(!is.na(allGO)) == 0) {
        return (NA)
        }
    if(!is.null(dropCodes)) {
        # if dropcodes exist, drop them.
        evidence <- sapply(allGO, function(x) x$Evidence)
        drop <- evidence %in% dropCodes
        allGO <- allGO[!drop]
        }
    # grab category
    category <- sapply(allGO, function(x) x$Ontology)
    #
    allGO <- allGO[category %in% ontology]

    if(length(allGO)==0) return (NA)
    return (unlist(unique(names(allGO)))) #return the GOIDs
}

# get gene ontology structure dataset.
data(go)
# calculate the information content for each go
IC <- descendants_IC(go)

# Re-used from ppiPre package (and gosim)
# gets all ancestors
getAncestors <- function(ont) {
    Ancestors <- switch(ont,
                        MF = "GOMFANCESTOR",
                        BP = "GOBPANCESTOR",
                        CC = "GOCCANCESTOR"
                        )

    return (eval(parse(text=Ancestors)))
}

# Re-used from ppiPre package (GetLatestCommonAncestor)
ComAnc<-function(GOID1,GOID2,ont, organism){
    rootCount <- max(IC[IC != Inf])
    IC["all"] = 0
    p1 <- IC[GOID1]/rootCount
    p2 <- IC[GOID2]/rootCount
    if(is.na(p1) || is.na(p2)) return (NA)
    if (p1 == 0 || p2 == 0) return (NA)

    ancestor1 <- unlist(getAncestors(ont)[[GOID1]])
    ancestor2 <- unlist(getAncestors(ont)[[GOID2]])
    # unlik ppipre the two GOIDs are always intersected (always the assumption that GOIDs will be different?)
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

# This will store simliarity pairs of GO terms to avoid repeated
# calculations, in particular to estimate similarity for a list of genes.
All_info_GO_Pairs <- data.frame(GO_1=0, GO_2=0, S_1=0)

# TopoICSim algorithm to calculate similarity between two GO terms.
TopoICSim_ti_tj <- function(GOID1, GOID2, ont, organism, WDG, GA, root){

  D_ti_tj_x <- c()
  COMANC<-ComAnc(GOID1, GOID2, ont, organism)
  if (length(COMANC)!=0 && !is.na(COMANC)){
    for (x in COMANC){
        # To identify all disjunctive common ancestors
        ImmadiateChildren_x <- switch(ont, MF = GOMFCHILDREN[[x]],
                                           BP = GOBPCHILDREN[[x]],
                                           CC = GOCCCHILDREN[[x]])
	if(x!="all"  &
           x!=root   &
           !is.na(x) &
           length(intersect(ImmadiateChildren_x, COMANC))==0){

           #Subgraph from two GO terms GOID1 and GOID2
            sg1 <- subGraph(c(get(GOID1, GA), GOID1), WDG)
            sg2 <- subGraph(c(get(GOID2, GA), GOID2), WDG)
            # Subgraph from a disjunctive common ancestor to root
	    sglca <- subGraph(c(get(x, GA), x), WDG)
            sglca <- WSG(sglca)
            wLP_x_root <- Longest_Path(sglca, x, root)

	    sg1 <- WSG(sg1)
	    sg1 <- igraph.to.graphNEL(sg1)
	    sp1 <- sp.between(sg1, GOID1, x)
	    ic_sp1 <- sp1[[1]]$length
	    length_sp1 <- length(sp1[[1]]$path_detail)

	    sg2 <- WSG(sg2)
	    sg2 <- igraph.to.graphNEL(sg2)
	    sp2 <- sp.between(sg2, GOID2, x)
	    ic_sp2 <- sp2[[1]]$length
	    length_sp2 <- length(sp2[[1]]$path_detail)

            IC_GOID_1 <- IC[GOID1][[1]]
	    IC_GOID_2 <- IC[GOID2][[1]]
            if (!is.na(IC_GOID_1) & IC_GOID_1!=0) ic_sp1 <- ic_sp1+(1/(2*IC_GOID_1))
            if (!is.na(IC_GOID_2) & IC_GOID_2!=0) ic_sp2 <- ic_sp2+(1/(2*IC_GOID_2))

	    wSP_ti_tj_x <- (ic_sp1+ic_sp2)*(length_sp1+length_sp2-2)
	    D_ti_tj_x <- c(D_ti_tj_x, wSP_ti_tj_x/wLP_x_root)
	    }
    }
  }
  if(!is.null(D_ti_tj_x)) {
	sim<-round(1-(atan(min(D_ti_tj_x,na.rm=TRUE))/(pi/2)), 3)
	All_info_GO_Pairs <- rbind(All_info_GO_Pairs, c(GOID1, GOID2, sim))
	assign("All_info_GO_Pairs", All_info_GO_Pairs,.GlobalEnv)
	return(sim)
	}
  else {
	All_info_GO_Pairs <- rbind(All_info_GO_Pairs, c(GOID1, GOID2, 0))
	assign("All_info_GO_Pairs", All_info_GO_Pairs,.GlobalEnv)
	return(0)
	}
}

# TopoICSim algorithm to calculate similarity between two genes or gene products.
TopoICSim <- function(gene1, gene2, ont="MF", organism="yeast", drop=NULL){
  # set ontology and organism
  ont <- match.arg(ont, c("MF", "BP", "CC"))
  organism <- match.arg(organism, c("human", "fly", "mouse",
                                    "rat", "yeast", "zebrafish",
                                    "worm", "arabidopsis", "ecolik12",
                                    "bovine", "canine", "anopheles",
                                    "ecsakai", "chicken", "chimp",
                                    "malaria", "rhesus", "pig", "xenopus"))
  # Initial definitions and sets based on organism and ontology type

  xx <- switch(ont, MF = toTable(GOMFPARENTS),
                    BP = toTable(GOBPPARENTS),
                    CC = toTable(GOCCPARENTS))

  GA <- switch(ont, MF = GOMFANCESTOR,
                    BP = GOBPANCESTOR,
                    CC = GOCCANCESTOR)

  root <- switch(ont, MF = "GO:0003674",
                      BP = "GO:0008150",
                      CC = "GO:0005575")
  WDG = ftM2graphNEL(as.matrix(xx[, 1:2]))

  go1 <- GetOntology(gene1, organism= organism, ontology= ont, dropCodes=drop)
  go2 <- GetOntology(gene2, organism= organism, ontology= ont, dropCodes=drop)

  if (sum(!is.na(go1)) == 0 || sum(!is.na(go2)) == 0) {
    return (list(geneSim=NA, GO1=go1, GO2=go2))
  }
  m <- length(go1)
  n <- length(go2)
  scores <- matrix(nrow=m, ncol=n)
  rownames(scores) <- go1
  colnames(scores) <- go2
  for( i in 1:m ) {
    for (j in 1:n) {
      if (go1[i]==go2[j]) scores[i,j] <- 1
      else {
      	DF<-subset(All_info_GO_Pairs, (All_info_GO_Pairs[,1] == go1[i] &
                                             All_info_GO_Pairs[,2] == go2[j]))
      	if (nrow(DF)==0) DF<-subset(All_info_GO_Pairs, (All_info_GO_Pairs[,1] == go2[j] &
                                                                    All_info_GO_Pairs[,2] == go1[i]))
      	if (nrow(DF)!=0) scores[i,j]<-as.numeric(DF[1,3])
      	else scores[i,j] <- TopoICSim_ti_tj(go1[i], go2[j], ont, organism, WDG, GA, root)
      	}
    }
  }
  if (!sum(!is.na(scores))) return(list(geneSim=NA, GO1=go1, GO2=go2))
  scores<-sqrt(scores)
  sim <- max(sum(sapply(1:m, function(x) {max(scores[x,], na.rm=TRUE)}))/m ,
             sum(sapply(1:n, function(x) {max(scores[,x], na.rm=TRUE)}))/n)
  sim <- round(sim, digits=3)
  return (list(geneSim=sim, GO1=go1, GO2=go2))
}

cat(rep("\n", 6))

# EXAMPLES

# Example 1
# Calculate TopoICSim between two human genes "219" and "8659"
start.time <- Sys.time()

TopoICSimValue<-TopoICSim("218", "501", ont="MF", organism="human", drop=NULL)
cat("############################     EXAMPLE 1     ############################")
cat("\n")
cat("\n")
cat("TopoICSim between two human genes 218 and 501 is: ")
cat("\n")
print(TopoICSimValue)

cat(rep("\n", 3))

# Example 2
# Calculate IntraSet similarity for first Pfam clan gene set
cat("############################     EXAMPLE 2     ############################")
cat("\n")
list1=c("126133","221","218","216","8854","220","219","160428","224","222","8659","501","64577","223","217","4329","10840","7915")
# Function to calculate IntraSet similarity
IntraSetSim <- function(List_Genes){
	IntraSim <- 0
        l <- length(List_Genes)
	for(i in 1:l){
		Gene1 <- toString(List_Genes[i])
		for(j in 1:l){
			Gene2 <- toString(List_Genes[j])
			TopoICSim_<-TopoICSim(Gene1, Gene2, ont="MF", organism="human", drop=NULL)[[1]]
			if(!is.na(TopoICSim_)) IntraSim <- IntraSim + TopoICSim_
			}
		}
	return(IntraSim /(l*l))
}

cat("IntraSetSim(CL0099.10)  =  ", IntraSetSim(list1), rep("\n", 6))

end.time <- Sys.time()
time.taken <- end.time - start.time

time.taken

