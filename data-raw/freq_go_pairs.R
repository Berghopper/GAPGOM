# In this file internal procedures will be defined for generating data, this is for "developers only".

# generate score matrix between most frequent occuring gos in all organisms/ontologies
.gen_semscor_topo_matrices <- function(old_listy = NULL, amount=2) {
  # old_listy is deprecated
  options(digits=22, stringsAsFactors = F)
  ids <- c("ENTREZ", "ENSEMBL")
  organisms <- c("mouse")#, "human")
  ontologies <- c("MF")#, "CC", "BP")
  combos <- expand.grid(ids, ontologies, organisms, stringsAsFactors = F)
  #colnames(combos) <- c("id", "organism", "ontology")
  res <- lapply(seq_len(nrow(combos)), function(i) {
    row <- combos[i,]
    id <- row[[1]]
    ontology <- row[[2]]
    organism <- row[[3]]
    tmp_godat <- GAPGOM::set_go_data(organism, ontology, computeIC = T)
    geneanno <- data.table(tmp_godat@geneAnno)
    counted <- geneanno[unique(geneanno[[1]]), .N, by=GO, on=.(ENTREZID)]
    top_gos <- counted[rev(order(counted$N)),][1:amount,][[1]]
    # now that we have the top gos, generate the matrix.
    custom <- setNames(as.list(top_gos), seq(top_gos))
    return(GAPGOM::topo_ic_sim_genes(ontology, organism, c(), c(), 
                                     custom_genes1 = custom, 
                                     custom_genes2 = custom, 
                                     go_data = tmp_godat, 
                                     use_precalculation = F, 
                                     garbage_collection = T, idtype = id)$AllGoPairs)
  })
  names(res) <- apply(combos, 1, function(row) {paste0(row, collapse = "_")})
  return(res)
}
res <- .gen_semscor_topo_matrices()
