#' GAPGOM internal - set_go_data()
#' 
#' Set GO data (this function purely makes choosing Bioconductor datasets a 
#' little easier)
#' 
#' @section Notes:
#' This function is an internal function and should not be called by the user.
#' @importFrom GOSemSim godata
set_go_data <- compiler::cmpfun(function(organism, ontology) {
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
  return(godata(species, ont = ontology, computeIC = TRUE))
})

#' strsplit(ft5_ranncom$table[!(is.na(ft5_ranncom$table$`entrezgene (genes) id associated with the transcript`) | ft5_ranncom$table$`entrezgene (genes) id associated with the transcript`==""),][1,5][[1]], "\\:")[[1]][2]

#' fantom5 data load.
#' @importFrom data.table fread
#' @examples 
#' ft5_r <- fantom_load_raw("/media/casper/USB_ccpeters/internship_thesis/data/hg19.cage_peak_phase1and2combined_counts.osc.txt2", verbose=T)
fantom_load_raw <- function(filepath, verbose = F) {
  old <- options(stringsAsFactors = F)
  on.exit(options(old), add = TRUE)
  # first parse fantom column variables.
  if (verbose) print("reading header variables...")
  con <- file(filepath, "r")
  translation_df <- data.frame(c(0,0,0))
  linecount <- 0
  while ( TRUE ) {
    line <- readLines(con, n = 1)
    if (!startsWith(line, "#")) {
      break
    } else {
      filtered_line <- sapply(strsplit(line, "##|\\[|\\]|\\="), 
                              function(x){x[!x ==""]})
      translation_df <- data.frame(translation_df, filtered_line)
    }
    linecount <- linecount + 1
  }
  translation_df <- as.data.frame(t(translation_df))
  rownames(translation_df) <- 1:nrow(translation_df)
  colnames(translation_df) <- c("type","key","value")
  close(con)
  if (verbose) print("DONE")
  if (verbose) print("loading table...")
  # load table
  fan <- fread(filepath, skip=linecount, showProgress = verbose)
  if (verbose) print("DONE")
  if (verbose) print("formatting column names...")
  # and touch up column names
  colnames(fan) <- sapply(colnames(fan), function(x){translation_df[translation_df$key==x,]$value})
  if (verbose) print("DONE")
  # return table+leftover metadata (header variables).
  return(list(table=fan, meta=translation_df[translation_df$type!="ColumnVariables" | translation_df$type!=0,]))
}

#' convert ft5 to expset
#' @importFrom Biobase ExpressionSet
fantom_to_expset <- function(fanraw, colselector, verbose = F) {
  fan <- fanraw$table
  meta <- fanraw$meta
  # try see if columnnames include any type of entrez id
  colnames_fan <- colnames(fan)
  exp <- regexec(".*entrez.*", colnames_fan)
  regex_result <- unlist(regmatches(colnames_fan, exp))
  if (length(regex_result)==0) {
    stop("ERROR: entrez id column not found!")
  } else if (length(regex_result) > 1) {
    print(length(regex_result))
    print(regex_result[1:100])
    stop("ERROR: mutliple entrez ids?")
  } else if (length(regex_result) == 1) {
    # if regex result contains only 1 string, move on.
    if (verbose) print("filtering out empty entrez ids...")
    fan <- fan[!(is.na(fan[[regex_result]]) | fan[[regex_result]]==""),]
    if (verbose) print("DONE")
    if (verbose) print("converting to expressionset...")
    expression_matrix <- as.matrix(fan[,colselector])
    rownames(expression_matrix) <- fan[[regex_result]]
    return(ExpressionSet(expression_matrix, ...=meta))
  }
}


#' convert into expressionset
fantom_load_expressionset <- function(filepath, colselector, verbose = F) {
  fanraw <- fantom_load_raw(filepath, verbose = verbose)
  return(fantom_to_expset(fanraw, colselector, verbose = verbose))
}


