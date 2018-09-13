
# R code to measure similarity two genes based on TopoICSim algorithm.
# Ehsani R, Drablos F. TopoICSim: A new semantic similarity measure based on gene ontology. In revision.
# Written by Rezvan Ehsani, <rezvanehsani74@gmail.com>.

###################
#rm(list = ls())

# load libs
library("GO.db")
library("graph")
library("RBGL")
library("igraph")
library("Rgraphviz")
library("GOSemSim")

options(warn=-1)
options(stringsAsFactors = FALSE)

# Set GO data
Set_GO_Data <-  function(organism, ontology) {
    species <- switch(organism,
                    human       = "org.Hs.eg.db",
                    fly         = "org.Dm.eg.db",
                    mouse       = "org.Mm.eg.db",
                    rat         = "org.Rn.eg.db",
                    yeast       = "org.Sc.sgd.db",
                    zebrafish   = "org.Dr.eg.db",
                    worm        = "org.Ce.eg.db",
                    arabidopsis = "org.At.tair.db",
                    ecolik12    = "org.EcK12.eg.db",
                    bovine	= "org.Bt.eg.db",
                    canine	= "org.Cf.eg.db",
                    anopheles	= "org.Ag.eg.db",
                    ecsakai	= "org.EcSakai.eg.db",
                    chicken	= "org.Gg.eg.db",
                    chimp	= "org.Pt.eg.db",
                    malaria	= "org.Pf.plasmo.db",
                    rhesus	= "org.Mmu.eg.db",
                    pig	        = "org.Ss.eg.db",
                    xenopus	= "org.Xl.eg.db"
    )
    godata(species, ont=ontology, computeIC=TRUE)
    }

# Common Ancestor for two GOIDs
ComAnc<-function(GOID1,GOID2,ont, organism, GA){
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

# This will store simliarity pairs of GO terms to avoid repeated
# calculations, in particular to estimate similarity for a list of genes.
All_info_GO_Pairs <- data.frame(GO_1=0, GO_2=0, S_1=0)

# TopoICSim algorithm to calculate similarity between two GO terms.
TopoICSim_ti_tj <- function(GOID1, GOID2, ont, organism, WDG, GA, root){

  D_ti_tj_x <- c()
  COMANC<-ComAnc(GOID1, GOID2, ont, organism, GA)
  if (length(COMANC)!=0 && !is.na(COMANC)){
    for (x in COMANC){
      # To identify all disjunctive common ancestors
      ImmadiateChildren_x <- switch(ont, MF = GOMFCHILDREN[[x]],
                                         BP = GOBPCHILDREN[[x]],
                                         CC = GOCCCHILDREN[[x]])
    	if(x!="all"  &
         x!=root   &
         !is.na(x) &
         length(intersect(ImmadiateChildren_x, COMANC))==0)
    	{
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
	# assign dataframe to global variable namespace
	assign("All_info_GO_Pairs", All_info_GO_Pairs,.GlobalEnv)
	return(0)
	}
}

# TopoICSim algorithm to calculate similarity between two genes or gene products.
TopoICSim_g1g2 <- function(gene1, gene2, ont="MF", organism="yeast", drop=NULL){
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
                      CC = "GO:0005575 ")
  WDG = ftM2graphNEL(as.matrix(xx[, 1:2]))

  goAnno <- GO_Data@geneAnno
  goAnno <- goAnno[!goAnno$EVIDENCE %in% drop,]
  go1 <- as.character(unique(goAnno[goAnno[,1] == gene1, "GO"]))
  go2 <- as.character(unique(goAnno[goAnno[,1] == gene2, "GO"]))

  if (sum(!is.na(go1)) == 0 || sum(!is.na(go2)) == 0) {
    return (list(geneSim=NA, GO1=go1, GO2=go2))
  }
  m <- length(go1)
  n <- length(go2)
  scores <- matrix(nrow=m, ncol=n)
  rownames(scores) <- go1
  colnames(scores) <- go2
  if(gene1==gene2) return(1)
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
  scores<-sqrt(scores) # Some scores seem unnecesary to calculate? Like, when they have the same GOID the similarity will always be 1.
			## YES, I ADDED A NEW LINE FOR THIS CASE ABOVE  //if(gene1==gene2) return(1)//
  #View(scores)
  sim <- max(sum(sapply(1:m, function(x) {max(scores[x,], na.rm=TRUE)}))/m ,
             sum(sapply(1:n, function(x) {max(scores[,x], na.rm=TRUE)}))/n)
  sim <- round(sim, digits=3)
  return(sim)
  #return (list(geneSim=sim, GO1=go1, GO2=go2))
}


# Calculate TopoICSim between two gene lists
TopoICSim <- function(GeneList1, GeneList2, ont="MF", organism="human", drop=NULL){
        GO_Data <- Set_GO_Data(organism=organism, ontology=ont)
        IC <- GO_Data@IC
        assign("GO_Data", GO_Data,.GlobalEnv)
        assign("IC", IC,.GlobalEnv)
        r1 <- length(GeneList1)
        r2 <- length(GeneList2)
        ScoreMatrix <- matrix(nrow=r1, ncol=r2)
        rownames(ScoreMatrix) <- GeneList1
        colnames(ScoreMatrix) <- GeneList2
        for(i in 1:r1){
            for(j in 1:r2){
                ScoreMatrix[i,j] <- TopoICSim_g1g2(GeneList1[i], GeneList2[j], ont, organism, drop)
                }
            }
        return(ScoreMatrix)
        }


cat(rep("\n", 6))

# EXAMPLES

# Example 1
# Calculate TopoICSim between two human genes "219" and "501"
TopoICSimValue<-TopoICSim("218", "501", ont="MF", organism="human", drop=NULL)
print("############################     EXAMPLE 1     ############################")
cat("\n")
cat("\n")
print("TopoICSim between two human genes 218 and 501 is: ")
cat("\n")
print(TopoICSimValue)

cat(rep("\n", 3))


# Example 2
# Calculate IntraSet similarity for first Pfam clan gene set
print("############################     EXAMPLE 2     ############################")
cat("\n")
list1=c("126133","221","218","216","8854","220","219","160428","224","222","8659","501","64577","223","217","4329","10840","7915")
# Function to calculate IntraSet similarity, Interset not implemented however?
## IN FACT HERE WE JUST PRESENTED TWO SIMPLE EXAMPLES THAT USER SEE HOW USE OF THE TopoICSim. WE WERE NOT GOING TO REPRODUCE ALL CALCULATIONS IN THE PAPER HERE.
IntraSetSim <- TopoICSim(list1, list1, ont="MF", organism="human", drop=NULL)
IntraSetSim
cat("IntraSetSim(CL0099.10)  =  ", mean(IntraSetSim))

cat(rep("\n", 6))

end.time <- Sys.time()
time.taken <- end.time - start.time

time.taken

