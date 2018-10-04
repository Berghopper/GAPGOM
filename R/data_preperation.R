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

#' GAPGOM - fantom_load_raw()
#'
#' Loads raw fantom5 data from file.
#'
#' This function loads raw fantom5 data and returns the resulting data.table/
#' data.frame object.
#'
#' @param filepath filename of fantom5 file.
#' @param verbose Switch to TRUE for extra messages. Default=FALSE
#'
#' @return The resulting datatable containing raw fantom5 data. (Most of the
#' time very large!)
#' 
#' @examples 
#' ft5 <- fantom_load_raw("./hg19.cage_peak_phase1and2combined_counts.osc.txt", 
#' verbose=T)
#' @importFrom data.table fread
#' @export
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
  if (verbose) print("loading dataframe...")
  # load df
  fan <- fread(filepath, skip=linecount, showProgress = verbose)
  if (verbose) print("DONE")
  if (verbose) print("formatting column names...")
  # and touch up column names
  colnames(fan) <- sapply(colnames(fan), function(x){translation_df[translation_df$key==x,]$value})
  if (verbose) print("DONE")
  # return df+leftover metadata (header variables).
  return(list(df=fan, meta=translation_df[translation_df$type!="ColumnVariables" | translation_df$type!=0,]))
}

#' GAPGOM - fantom_to_expset()
#'
#' Convert raw data.table/data.frame fantom5 object to a proper ExpressionSet.
#'
#' This function converts fantom5 data and converts it into an ExpresionSet.
#' This ExpressionSet is then returned.
#' This function only accepts the RLE normalized data!
#'
#' @param fanraw raw data.table object from the fantom_load_raw() function.
#' @param filter Filter, this causes only entries to be added that have an
#' entrez ID. Normally this should be left on default (TRUE) because all
#' algorithms in this library need the entrez IDs for translation. 
#' @param verbose Switch to TRUE for extra messages. Default=FALSE
#'
#' @return The resulting ExpressionSet contains the original data. The 
#' expressiondata can be found under ExpressionSet@assayData$exprs
#' Other information (first 6 info columns) can be found under;
#' ExpressionSet@featureData@data
#' 
#' @examples 
#' ft5 <- fantom_load_raw("./mm9.cage_peak_phase1and2combined_tpm_ann.osc.txt", 
#' verbose = T)
#' expset <- fantom_to_expset(ft5, verbose = T)
#' @importFrom Biobase ExpressionSet annotatedDataFrameFrom
#' @export
fantom_to_expset <- function(fanraw, filter = T, verbose = F) {
  fan <- fanraw$df
  meta <- fanraw$meta
  if (filter) {
    if (verbose) print("filtering out empty entrez ids...")
    fan <- fantom_filter_entrez(fan)
  } 
  if (verbose) print("DONE")
  if (verbose) print("converting to expressionset...")
  expression_matrix <- as.matrix(fan[,7:ncol(fan)])
  rownames(expression_matrix) <- fan$`CAGE peak id`
  featuredat <- as.data.frame(fan[,2:6])
  rownames(featuredat) <- fan$`CAGE peak id`
  return(ExpressionSet(expression_matrix, 
                       featureData = new("AnnotatedDataFrame", 
                                         data=featuredat),...=meta)) # #6 issue
}

#' GAPGOM internal - fantom_filter_entrez()
#'
#' filter on entrez id with raw fantom5 matrix.
#' 
#' @section Notes:
#' This function is an internal function and should not be called by the user.
fantom_filter_entrez <- function(fan) {
  return(fan[!(
    is.na(fan$`entrezgene (genes) id associated with the transcript`) | 
      fan$`entrezgene (genes) id associated with the transcript`==""),])
}

#' GAPGOM - fantom_load_expressionset()
#'
#' Loads and converts fantom5 data to an ExpressionSet
#'
#' This function loads fantom5 data and then converts it into an ExpresionSet.
#' This ExpressionSet is then returned. When you don't want to instantly load
#' and convert, use load_fantom_raw() instead and then convert it with
#' fantom_to_expset(). The standard annotation columns always need to be kept!
#' (First 6 columns). This is also true for at least 1 or more expression
#' column(s). This function only accepts the RLE normalized data!
#'
#' @param filepath filename of fantom5 file.
#' @param filter Filter, this causes only entries to be added that have an
#' entrez ID. Normally this should be left on default (TRUE) because all
#' algorithms in this library need the entrez IDs for translation. 
#' @param verbose Switch to TRUE for extra messages. Default=FALSE
#'
#' @return The resulting ExpressionSet contains the original data. The 
#' expressiondata can be found under ExpressionSet@assayData$exprs
#' Other information (first 6 info columns) can be found under;
#'  ExpressionSet@featureData@data
#' 
#' @examples
#' fantom_file <- "./mm9.cage_peak_phase1and2combined_tpm_ann.osc.txt"
#' expset <- fantom_load_expression_set(fantom_file)
#' View(expset)
#' @export
fantom_load_expressionset <- function(filepath, filter = T, verbose = F) {
  fanraw <- fantom_load_raw(filepath, verbose = verbose)
  return(fantom_to_expset(fanraw, filter = filter, verbose = verbose))
}

#' GAPGOM - fantom_download()
#'
#' Downloads and unpacks fantom5 data of either human or mouse.
#'
#' This function downloads the whole fantom5 dataset and unpacks it. automatic
#' unpacking can be turned off
#'
#' @param down_dir The folder/directory to download the fantom5 file to.
#' @param organism Either "mouse" or "human". FANTOM5 only has these two as
#' fully annotated+tpm normalized datasets. 
#' @param unpack Default=TRUE, if set to TRUE, file will be unpacked 
#' automatically.
#' @param noprompt Default=FALSE, user prompt, if you are 100% sure you want to
#' download and have enough space, set this to TRUE.
#'
#' @return The resulting filename/location of the file or NULL if cancelled.
#' 
#' @examples
#' fantom_file <- fantom_download("./", organism = "mouse", noprompt = T)
#' @importFrom GEOquery gunzip
#' @export
fantom_download <- function(down_dir, organism="human", unpack = T, 
                            noprompt = F) {
  baseurl <- "http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/"
  url <- switch(organism, 
                "human" = paste0(
                  baseurl,
                  "hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz"
                ),
                "mouse" = paste0(
                  baseurl,
                  "mm9.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz"
                ), 
                stop(paste0("INCORRECT ORGANISM; \"", organism, "\" 
                                  SELECTED! --> INVALID OPTION"))
  )
  filename <- strsplit(url, "/")[[1]][length(strsplit(url, "/")[[1]])]
  full_filename <- paste0(down_dir,"/",filename)
  # check for user prompt if selected
  if(!noprompt) {
    rawsize <- switch(organism, 
                      "human"=c(792, 2200), 
                      "mouse"=c(437, 1200), 
                      stop(paste0("INCORRECT ORGANISM; \"", organism, "\" 
                                  SELECTED! --> INVALID OPTION")))
    size_selector <- 1
    if (unpack) {
      size_selector <- 2
    }
    size <- rawsize[size_selector]
    cat(paste0(
      "Are you sure you want to download the \"", 
      filename, "\" fantom5 file?\nFINAL SIZE; ",size,"MB\n",
      "Answer \"yes\" or \"y\": "))
    answer <- tolower(readline())
    if (answer == "y"|answer == "yes") {
      print("Starting download!")
    } else {
      print("Download cancelled!")
      # return null to cancel
      return(NULL)
    }
  } 
  download.file(url, full_filename, "auto")
  if (unpack) {
    gunzip(full_filename, overwrite = T)
    full_filename <- substr(full_filename, 1, len(full_filename)-3)
  }
  return(full_filename)
}
