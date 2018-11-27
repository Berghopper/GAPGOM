# File containing user argument checks (SHOULD BE EXPANDED)

### OTHER CHECKS / GLOBAL CHECKS

.check_ifclass <- function(obj, 
               classname, 
               objname, 
               match_case = FALSE, 
               accept_null = TRUE) {
  obj_classname <- class(obj)[1]
  if (!match_case) {
    obj_classname <- tolower(obj_classname)
    classname <- tolower(classname)
  }
  
  if(obj_classname!=classname) {
    if(is.null(obj) && accept_null) {
      return(TRUE)
    }
    message(paste0("\"", objname, "\" should be a \"", classname, "\"!"))
    return(FALSE)
  }
  return(TRUE)
}

.check_gos <- function(gos) {
  if(class(gos)!="character") {
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
  message("One of the custom genes is not a list, is missing names or has 
          incorrect gos!")
  return(FALSE)
}

.check_idtype <- function(idtype, organism) {
  species <- .organism_to_species_lib(organism)
  if (!is.null(species)) {
    if (idtype %in% keytypes(eval(parse(text=species)))) {
      return(TRUE)
    }
  }
  message(paste0("idtype is incorrect! Pick one of the following; ", 
                 paste0(keytypes(eval(parse(text=species))), collapse = ", ")))
  return(FALSE)
}

.check_organism <- function(organism) {
  if(is.null(.organism_to_species_lib(organism))) {
    message("Incorrectly specified organism!")
    return(FALSE)
  }
  return(TRUE)
}

.check_ontology <- function(ontology) {
  ontologies <- c("MF", "BP", "CC")
  if(toupper(ontology) %in% ontologies) {
    return(TRUE)
  }
  message(paste0("Incorrectly specific ontology! Must be one of the 
                 following: ",
                 message(paste0("\"", 
                                paste0(ontologies, collapse="\", \""), 
                                "\""))))
  return(FALSE)
}

.check_method <- function(method) {
  methods <- c("pearson", "spearman", "kendall", "fisher", "sobolev", "combine")
  if(tolower(method) %in% methods) {
    return(TRUE)
  }
  message(paste0("Incorrectly specific method! Must be one of the following: ",
                 message(paste0("\"", 
                                paste0(methods, collapse="\", \""), 
                                "\""))))
  return(FALSE)
}