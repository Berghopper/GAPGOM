# File containing user argument checks (SHOULD BE EXPANDED)

#' GAPGOM internal - .topo_ic_sim_argchecks_genes()
#' 
#' This function is an internal function and should not be called by the user.
#' 
#' This function checks arguments for the topoicsim gene variant. This function
#' is temporary and will be replaced with other checks later on.
#' 
#' @keywords internal
.topo_ic_sim_argcheck_genes <- compiler::cmpfun(function(ontology, organism, genes1, genes2) {
  # set ontology and organism
  ontology <- match.arg(ontology, c("MF", "BP", "CC"))
  organism <- match.arg(organism, c("human",
                                    "fly",
                                    "mouse",
                                    "rat",
                                    "yeast",
                                    "zebrafish",
                                    "worm",
                                    "arabidopsis",
                                    "ecolik12",
                                    "bovine",
                                    "canine",
                                    "anopheles",
                                    "ecsakai",
                                    "chicken",
                                    "chimp",
                                    "malaria",
                                    "rhesus",
                                    "pig",
                                    "xenopus"))
  if (length(genes1) < 1 || length(genes2) < 1) {
    stop("Not enough genes specified!")
  }
})