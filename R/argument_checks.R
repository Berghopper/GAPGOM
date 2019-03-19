# File containing user argument checks (SHOULD BE EXPANDED)

### OTHER CHECKS / GLOBAL CHECKS

#' GAPGOM internal - Global argument checks 
#' 
#' These functions are internal functions and should not be called by the user.
#' 
#' These check make sure that user input is correct for both TopoICSim and 
#' lncRNApred. They mainly include very generalistic checks.
#'  
#' @section Notes:
#' Internal functions used in both algorithms.
#'
#' @return output is different on a case-to-case basis, but is a boolean 
#' determining if the input is correct.
#'
#' @name glb_checks
#' @keywords internal
NULL

#' Checks if object is the given classtype. prints message if incorrect.
#' @rdname glb_checks
#' @importFrom methods is
.check_ifclass <- function(obj, classname, objname, accept_null = TRUE) {
  if(!is(obj, classname)) {
    if(is.null(obj) & accept_null) {
      return(TRUE)
    }
    message("\"", objname, "\" should be a \"", classname, "\"!")
    return(FALSE)
  }
  return(TRUE)
}

#' Checks if given Go strings are correct.
#' @rdname glb_checks
#' @importFrom methods is
.check_gos <- function(gos) {
  if(!is(gos, "character")) {
    message("GOs should be character class!")
    return(FALSE)
  }
  regex_str <- "GO:\\d*$"
  go_matches <- unlist(regmatches(gos, regexec(regex_str, gos)), FALSE, FALSE)
  if(length(gos) == length(go_matches)) {
    return(TRUE)
  }
  message("One or multiple GOs are not correct!")
  return(FALSE)
}

#' Checks if a custom gene is correctly defined. prints message if incorrect.
#' @rdname glb_checks
.check_custom_gene <- function(cus_gen) {
  if(.check_ifclass(cus_gen, "list", "custom_gene")) {
    if(!is.null(names(cus_gen))) {
      if(.check_gos(unlist(cus_gen, FALSE, FALSE))) {
        return(TRUE)  
      }
    }
  }
  if(is.null(cus_gen)) {
    return(TRUE)
  }
  message("One of the custom genes is not a list, is missing names or has ",
    "incorrect gos!")
  return(FALSE)
}

#' Checks if idtype is correct and prints message if incorrect.
#' @rdname glb_checks
.check_idtype <- function(idtype, organism) {
  species <- .organism_to_species_lib(organism)
  if (!is.null(species)) {
    if (idtype %in% keytypes(eval(parse(text=species)))) {
      return(TRUE)
    }
  }
  message("idtype is incorrect! Pick one of the following: ", 
    paste0(keytypes(eval(parse(text=species))), collapse = ", "))
  return(FALSE)
}

#' Checks if organism can be parsed. prints message if incorrect.
#' @rdname glb_checks
.check_organism <- function(organism) {
  if(is.null(.organism_to_species_lib(organism))) {
    message("Incorrectly specified organism!")
    return(FALSE)
  }
  return(TRUE)
}

#' Checks if ontology can be parsed. prints message if incorrect.
#' @rdname glb_checks
.check_ontology <- function(ontology) {
  ontologies <- c("MF", "BP", "CC")
  if(toupper(ontology) %in% ontologies) {
    return(TRUE)
  }
  message("Incorrectly specific ontology! Must be one of the following: ",
    "\"", paste0(ontologies, collapse="\", \""), "\"")
  return(FALSE)
}

#' Checks if method can be parsed. prints message if incorrect.
#' @rdname glb_checks
.check_method <- function(method) {
  methods <- c("pearson", "spearman", "kendall", "fisher", "sobolev", "combine")
  if(tolower(method) %in% methods) {
    return(TRUE)
  }
  message("Incorrectly specific method! Must be one of the following: ",
    "\"", paste0(methods, collapse="\", \""), "\"")
  return(FALSE)
}
