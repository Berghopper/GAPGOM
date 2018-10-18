#' GAPGOM internal - set_go_data()
#' 
#' This function is an internal function and should not be called by the user.
#'
#' Set GO data (this function purely makes choosing Bioconductor datasets a 
#' little easier)
#'
#' @section Notes:
#' Internal function used in multiple functions of topoICSim.
#' 
#' @param organism where to be scanned genes reside in, this option
#' is neccesary to select the correct GO DAG. Options are based on the org.db
#' bioconductor package;
#' http://www.bioconductor.org/packages/release/BiocViews.html#___OrgDb
#' Following options are available: "fly", "mouse", "rat", "yeast",
#' "zebrafish", "worm", "arabidopsis", "ecolik12", "bovine", "canine",
#' "anopheles", "ecsakai", "chicken", "chimp", "malaria", "rhesus", "pig",
#' "xenopus". Fantom5 data only has "human" and "mouse" available depending
#' on the dataset.
#' @param ontology desired ontology to use for prediction. One of three;
#' "BP" (Biological process), "MF" (Molecular function) or "CC"
#' (Cellular Component). Cellular Component is not included with the package's
#' standard data and will thus yield no results.
#' @param copmuteIC whether to compute Information Content.
#' 
#' @return return godata as from GoSemSim
#' 
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

###TOPOICSIM FUNCTIONS

#' GAPGOM internal - .prepare_score_matrix_topoicsim()
#'
#' This function is an internal function and should not be called by the user.
#'
#' Prepared score dataframe for the semantic measures
#'
#' @section Notes:
#' Internal function used in multiple functions of topoICSim.
#'
#' @param vec1 vector1 with arbitrary lookup names
#' @param vec2 vector2 with arbitrary lookup names.
#' 
#' @return The score matrix with names and NA's.
#' @importFrom Matrix Matrix
.prepare_score_matrix_topoicsim <- compiler::cmpfun(function(vec1, vec2) {
  score_matrix <- Matrix(nrow = length(vec1), 
         ncol = length(vec2), 
         dimnames = list(
           vec1, vec2))
  score_matrix <- .set_identical_items(score_matrix)
  return(score_matrix)
})

#' GAPGOM internal - .set_identical_items()
#'
#' This function is an internal function and should not be called by the user.
#'
#' Sets identical items within a score matrix, the similarity will always be 1.
#'
#' @section Notes:
#' Internal function used in various topoclsim algorithms.
#'
#' @param score_matrix score matrix for topoclsim
#' 
#' @return Same matrix with correct sets set to 1.
.set_identical_items <- compiler::cmpfun(function(score_matrix) {
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
})

#' GAPGOM internal - .unique_combos()
#'
#' This function is an internal function and should not be called by the user.
#'
#' Generates unique combinations of two vectors
#'
#' @section Notes:
#' Internal function used in various topoclsim algorithms.
#'
#' @param v1 first vector to be combined with second.
#' @param v2 second vector to be combined with first.
#' 
#' @return datatable with unique combinations of the vectors excluding items
#' with same name and double combinations
#' 
#' @import data.table
#' @importFrom plyr .
.unique_combos <- compiler::cmpfun(function(v1, v2) {
  intermediate <- unique(CJ(v1, v2)[V1 > V2, c("V1", "V2") := .(V2, V1)])
  return(intermediate[V1 != V2])
})

### LNCRNAPRED FUNCTIONS

#' GAPGOM internal - .prepare_score_df()
#'
#' This function is an internal function and should not be called by the user.
#'
#' Prepared score dataframe for the enrichment analysis.
#'
#' @section Notes:
#' Internal function used in prediction_function().
#'
#' @param original_ids rowname ID's for selected expression data.
#' @param score score vector
#' @param gene_id Gene ID originally selected in the prediction function. 
#' (gets filtered)
#' 
#' @return The score dataframe with ensmbl ID's
#' 
#' @importFrom stats na.omit
.prepare_score_df <- compiler::cmpfun(function(original_ids, score, gene_id) {
  # score_df is the 'score' (correlation/metric) dataframe against target gene.
  # combine the dataframe if combine is selected, otherwise just add regular score.
  score_df <- data.frame(original_ids, score)
  score_df <- na.omit(score_df)
  # filter gene ID's (Select everything except the chosen gene).
  score_df <- score_df[(score_df[, 1] != gene_id), ]
  return(score_df)
})
