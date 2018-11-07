# In this file internal procedures will be defined for generating data, this is for "developers only".

# generate score matrix between most frequent occuring gos in all organisms/ontologies
.gen_semscor_topo_mat <- function(old_listy = NULL, amount=10) {
  library(GO.db)
  library(data.table)
  library(GOSemSim)
  library(graph)
  library(igraph)
  library(GAPGOM)
  library(pryr)
  options(digits=22)
  organisms <- c("mouse", "human")
  ontologies <- c("MF", "CC", "BP")
  all_scores <<- list()
  for (onto in ontologies) {
    for (organism in organisms) {
      .gen_single_semscore_topo_mat(organism=organism, onto=onto, all_scores=all_scores, amount=amount, old_df=old_listy[[paste0(onto,"_",organism)]])
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
  top_100_gos <- counted[rev(order(counted$N)),][1:amount,][[1]]
  scores <- GAPGOM:::.prepare_score_matrix_topoicsim(top_100_gos, top_100_gos, sparse=T)
  IC <-tmp_godat@IC
  unique_pairs <- GAPGOM:::.unique_combos(top_100_gos, top_100_gos)
  
  #existing_pairs <- NULL
  if (!is.null(old_df)) {
    incommon <- top_100_gos[top_100_gos %in% rownames(old_df)]
    existing_pairs <- GAPGOM:::.unique_combos(incommon, incommon)
    print(paste0("Detected old values! adding ", paste0(nrow(existing_pairs))," values..."))
    unique_pairs <- dplyr::anti_join(unique_pairs, existing_pairs , by=c("V1","V2"))
    old_dim_rows <- dim(old_df)[1]
    old_dim_cols <- dim(old_df)[2]
    scores[1:old_dim_rows, 1:old_dim_cols] <- old_df
  }
  print(paste0("Adding ", nrow(unique_pairs), " similarties."))
  print(Sys.time())
  pb <- txtProgressBar(min = 0, max = nrow(unique_pairs), style = 3)
  for (i in seq_len(nrow(unique_pairs))) {
    pair <- unique_pairs[i,]
    setTxtProgressBar(pb, i)
    if ((i %% 500) == 0) {
      gc()
    }
    go1 <- pair[[1]]
    go2 <- pair[[2]]
    scores <- GAPGOM:::.set_values(go1, go2, scores,
                                    GAPGOM:::.topo_ic_sim_titj(go1,
                                                               go2,
                                                               onto,
                                                               organism,
                                                               weighted_dag,
                                                               go_annotation,
                                                               root,
                                                               IC))
  }
  print(Sys.time())
  print(mem_used())
  print("Clearing memory...")
  gc()
  print(mem_used())
  all_scores[[paste0(onto,"_",organism)]] <<- scores
  save(all_scores, file = filename, compress = "xz", compression_level = 9)
}

# between each calculation, session is restarted to clear ram. R has a lot of
# ram leaks for some reason...
.gen_semscor_topo_mat(amount=150)
all_scores_bckp <- all_scores
.gen_semscor_topo_mat(amount=250, old_listy=all_scores_bckp)
all_scores_bckp <- all_scores
.gen_semscor_topo_mat(amount=300, old_listy=all_scores_bckp)
all_scores_bckp <- all_scores
.gen_semscor_topo_mat(amount=350, old_listy=all_scores_bckp)