# In this file internal procedures will be defined for generating data, this is for "developers only".

# generate score matrix between most frequent occuring gos in all organisms/ontologies
.gen_semscor_topo_mat <- function(old_listy = NULL, amount=10) {
  library(GO.db)
  library(data.table)
  library(GOSemSim)
  library(graph)
  library(igraph)
  library(GAPGOM)
  options(digits=22)
  organisms <- c("mouse", "human")
  ontologies <- c("MF", "CC", "BP")
  all_scores <<- list()
  for (onto in ontologies) {
    for (organism in organisms) {
      .gen_single_semscore_topo_mat(organism=organism, onto=onto, all_scores=all_scores, amount=amount)
    }  
  }
  
}

# generate a single score matrix and save the progress
.gen_single_semscore_topo_mat <- function(organism, onto, filename = "~/all.RData", old_df = NULL, amount = 10, all_scores = list()) {
  xx_parents <- switch(onto, MF = toTable(GOMFPARENTS),
                       BP = toTable(GOBPPARENTS), CC = toTable(GOCCPARENTS))
  
  go_annotation <- switch(onto, MF = GOMFANCESTOR, BP = GOBPANCESTOR,
                          CC = GOCCANCESTOR)
  
  root <- switch(onto, MF = "GO:0003674", BP = "GO:0008150",
                 CC = "GO:0005575")
  weighted_dag <- ftM2graphNEL(as.matrix(xx_parents[, 1:2]))
  tmp_godat <- GAPGOM:::.set_go_data(organism, onto, computeIC = T)
  geneanno <- data.table(tmp_godat@geneAnno)
  counted <- geneanno[unique(geneanno[[1]]), .N, by=GO, on=.(ENTREZID)]
  top_100_gos <<- counted[rev(order(counted$N)),][1:amount,][[1]]
  print(top_100_gos)
  scores <- GAPGOM:::.prepare_score_matrix_topoicsim(top_100_gos, top_100_gos)
  IC <-tmp_godat@IC
  unique_pairs <- GAPGOM:::.unique_combos(top_100_gos, top_100_gos)
  
  if (!is.null(old_df)) {
    existing_pairs <- GAPGOM:::.unique_combos(top_100_gos[top_100_gos %in% rownames(old_df)])
    dplyr::anti_join(by=c("V1","V2"))
    df1[!(df1$name %in% df2$name),]
  }
  # if (is.null(old_listy)) {
  #   old_scores <- NA
  # } else {
  #   old_scores <- old_listy[[paste0(onto,"_",organism)]]
  # }
  print(Sys.time())
  pb <- txtProgressBar(min = 0, max = nrow(unique_pairs), style = 3)
  proggy <- 0
  #print(nrow(unique_pairs))
  apply(unique_pairs, 1, function(pair) {
    setTxtProgressBar(pb, proggy)
    proggy <<- proggy + 1
    #print(proggy)
    if ((proggy %% 500) == 0) {
      #print("gc!")
      gc()
    }
    go1 <- pair[1]
    go2 <- pair[2]
    # if this is not the case (row is not present), then run topo_ic_sim 
    # between 2 terms.
    # if (go1 %in% rownames(old_scores) && go2 %in% rownames(old_scores)) {
    # scores <<- .set_values(go1, go2, scores, old_scores[go1, go2])
    # } else {
    scores <<- GAPGOM:::.set_values(go1, go2, scores,
                                    GAPGOM:::.topo_ic_sim_titj(go1,
                                                               go2,
                                                               onto,
                                                               organism,
                                                               weighted_dag,
                                                               go_annotation,
                                                               root,
                                                               IC))
    #gc()
    # }
    
  })
  print(Sys.time())
  gc()
  all_scores[[paste0(onto,"_",organism)]] <<- scores
  save(all_scores, file = filename, compress = "xz", compression_level = 9)
}

.gen_semscor_topo_mat(amount = 10)
