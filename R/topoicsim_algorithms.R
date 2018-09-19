# some very ugly global vars should be fixed;
# IC (present in "TopoICSim_ti_tj", "Longest_Path", "WSG", "ComAnc")
# All_info_GO_Pairs (present in "TopoICSim_g1g2", "TopoICSim_ti_tj")
# GO_Data (present in "TopoICSim_g1g2")

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
        # get the subset of the GO pair.
        DF<-subset(All_info_GO_Pairs, (All_info_GO_Pairs[,1] == go1[i] &
                                         All_info_GO_Pairs[,2] == go2[j]))
        # now check if they are present, if 0 rows are returned, then check the inverse combination
        if (nrow(DF)==0) DF<-subset(All_info_GO_Pairs, (All_info_GO_Pairs[,1] == go2[j] &
                                                          All_info_GO_Pairs[,2] == go1[i]))
        # if this is not returning 0, set score to whatever it already is in the scores df.
        if (nrow(DF)!=0) scores[i,j]<-as.numeric(DF[1,3])
        # if this is not the case (row is not present), then run TopoICSim between 2 terms.
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
