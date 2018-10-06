#' GAPGOM internal - set_go_data()
#' 
#' Set GO data (this function purely makes choosing Bioconductor datasets a 
#' little easier)
#' 
#' @section Notes:
#' This function is an internal function and should not be called by the user.
#' @importFrom GOSemSim godata
.set_go_data <- compiler::cmpfun(function(organism, ontology, computeIC = T) {
  species <- switch(organism, human = "org.Hs.eg.db",
                    fly = "org.Dm.eg.db",
                    mouse = "org.Mm.eg.db",
                    rat = "org.Rn.eg.db",
                    yeast = "org.Sc.sgd.db",
                    zebrafish = "org.Dr.eg.db",
                    worm = "org.Ce.eg.db",
                    arabidopsis = "org.At.tair.db",
                    ecolik12 = "org.EcK12.eg.db",
                    bovine = "org.Bt.eg.db",
                    canine = "org.Cf.eg.db",
                    anopheles = "org.Ag.eg.db",
                    ecsakai = "org.EcSakai.eg.db",
                    chicken = "org.Gg.eg.db",
                    chimp = "org.Pt.eg.db",
                    malaria = "org.Pf.plasmo.db",
                    rhesus = "org.Mmu.eg.db",
                    pig = "org.Ss.eg.db",
                    xenopus = "org.Xl.eg.db")
  return(godata(species, ont = ontology, computeIC = computeIC)) #KEYS SUPPORT!
})


.prepare_score_matrix_topoicsim <- function(vec1, vec2) {
  score_matrix <- matrix(nrow = length(vec1), 
         ncol = length(vec2), 
         dimnames = list(
           vec1, vec2))
  score_matrix <- .set_identical_items(score_matrix)
  return(score_matrix)
}


.set_identical_items <- function(score_matrix) {
  # set all matching names of matrix to 1 (Same genes).
  matched_by_row <- match(rownames(score_matrix), colnames(score_matrix))
  #print(matched_by_row)
  #matched_by_row <- matched_by_row[!is.na(matched_by_row)]
  lapply(which(!is.na(matched_by_row[seq_along(matched_by_row)])), function(i) {
      row <- i
      col <- matched_by_row[i]
      score_matrix[row, col] <<- 1.0
  })
  return(score_matrix)
}

#' @import data.table
#' @importFrom plyr .
.unique_combos <- function(v1, v2) {
  intermediate <- unique(CJ(v1, v2)[V1 > V2, c("V1", "V2") := .(V2, V1)])
  return(intermediate[V1 != V2])
}