#' GAPGOM - fantom_load_raw()
#'
#' Loads raw fantom5 data from file.
#'
#' This function loads raw fantom5 data and returns the resulting data.table/
#' data.frame object.
#'
#' @param filepath filename of fantom5 file.
#' @param verbose Switch to TRUE for extra messages. Default=FALSE
#' @param example Boolean switch for R CMD Check (NOT MEANT TO BE TURNED ON FOR
#' END-USERS).
#'
#' @return The resulting datatable containing raw fantom5 data. (Most of the
#' time very large!)
#' 
#' @examples 
#' fantom_file <- fantom_download(organism = "mouse", noprompt = TRUE)
#' ft5 <- fantom_load_raw(fantom_file, verbose = TRUE, example = TRUE)
#' @importFrom data.table fread
#' @export
fantom_load_raw <- function(filepath, verbose = FALSE, example = FALSE) {
  old <- options(stringsAsFactors = FALSE)
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
      filtered_line <- vapply(strsplit(line, "##|\\[|\\]|\\="), 
                              function(x){
                                x <- x[!x ==""]
                                x <- x[!x ==" "]
                                return(x)
                              }, character(3))
      translation_df <- data.frame(translation_df, filtered_line)
    }
    linecount <- linecount + 1
  }
  translation_df <- as.data.frame(t(translation_df))
  rownames(translation_df) <- seq_len(nrow(translation_df))
  colnames(translation_df) <- c("type","key","value")
  close(con)
  if (verbose) print("DONE")
  if (verbose) print("loading dataframe...")
  # load df
  if (example) {
    fan <- fread(filepath, skip=linecount, showProgress = verbose, nrows = 10)
  } else {
    fan <- fread(filepath, skip=linecount, showProgress = verbose)
  }
  if (verbose) print("DONE")
  if (verbose) print("formatting column names...")
  # and touch up column names
  colnames(fan) <- vapply(colnames(fan), function(x){
    translation_df[translation_df$key==x,]$value}, character(1))
  if (verbose) print("DONE")
  # return df+leftover metadata (header variables).
  return(list(df=fan, 
              meta=translation_df[translation_df$type!="ColumnVariables" | 
                                    translation_df$type!=0,]))
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
#' @param species either "human" or "mouse". This is important because both
#' datasets have different metadata/stats
#' @param filter Filter, this causes only entries to be added that have an
#' entrez ID. Normally this should be left on default (TRUE) because all
#' algorithms in this library need the entrez IDs for translation. 
#' @param verbose Switch to TRUE for extra messages. Default=FALSE
#'
#' @return The resulting ExpressionSet contains the original data. The 
#' expressiondata can be found under assayData(ExpressionSet)[["exprs"]]
#' Other information (first 6 info columns) can be found under;
#' pData(featureData(ExpressionSet))
#' 
#' @examples 
#' fantom_file <- fantom_download(organism = "mouse", noprompt = TRUE)
#' ft5 <- fantom_load_raw(fantom_file, verbose = TRUE, example = TRUE)
#' expset <- fantom_to_expset(ft5, "mouse", verbose = TRUE)
#' @importFrom Biobase ExpressionSet annotatedDataFrameFrom
#' @importFrom methods new
#' @export
fantom_to_expset <- function(fanraw, species, filter = TRUE, verbose = FALSE) {
  fan <- fanraw$df
  meta <- fanraw$meta
  if (filter) {
    if (verbose) print("filtering out empty entrez ids...")
    if (verbose) print("DONE")
    fan <- .fantom_filter_entrez(fan)
  } 
  icol <- switch (species,
          human=8,
          mouse=7)
  if (verbose) print("converting to expressionset...")
  expression_matrix <- as.matrix(fan[,icol:ncol(fan)])
  rownames(expression_matrix) <- fan$`CAGE peak id`
  featuredat <- as.data.frame(fan[,2:(icol-1)])
  rownames(featuredat) <- fan$`CAGE peak id`
  expset <- ExpressionSet(expression_matrix, 
                          featureData = new("AnnotatedDataFrame", 
                                            data=featuredat),...=meta) 
  # issue #6
  if (verbose) print("DONE")
  return(expset)
}

#' GAPGOM internal - .fantom_filter_entrez()
#'
#' filter on entrez id with raw fantom5 matrix and only keep rows with entrez 
#' ids.
#' 
#' @param fan fantom5 data.table/data.frame object
#' 
#' @section Notes:
#' This function is an internal function and should not be called by the user.
#' 
#' @return output is different on a case-to-case basis
#' 
#' @keywords internal
.fantom_filter_entrez <- function(fan) {
  return(fan[!(
    is.na(fan$`entrezgene (genes) id associated with the transcript`) | 
      fan$`entrezgene (genes) id associated with the transcript`==""),])
}

#' GAPGOM - fantom_download()
#'
#' Downloads and unpacks fantom5 data of either human or mouse.
#'
#' This function downloads the whole fantom5 dataset and unpacks it. automatic
#' unpacking can be turned off. The file is downloaded to a special caching
#' directory, which will be returned on exit.
#'
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
#' fantom_file <- fantom_download(organism = "mouse", noprompt = TRUE)
#' 
#' @importFrom GEOquery gunzip
#' @importFrom BiocFileCache BiocFileCache bfcrpath bfcquery bfcnew
#' @export
fantom_download <- function(organism="human", unpack = TRUE, noprompt = FALSE) {
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
                stop("INCORRECT ORGANISM; \"", organism, "\" ",
                  "SELECTED! --> INVALID OPTION")
  )
  filename <- strsplit(url, "/")[[1]][length(strsplit(url, "/")[[1]])]
  # check for user prompt if selected
  if(!noprompt) {
    rawsize <- switch(organism, 
                      "human"=c(792, 792+2200), 
                      "mouse"=c(437, 437+1200), 
                      stop("INCORRECT ORGANISM; \"", organism, "\" ",
                        "SELECTED! --> INVALID OPTION"))
    size_selector <- 1
    if (unpack) {
      size_selector <- 2
    }
    size <- rawsize[size_selector]
    cat(paste0(
      "Are you sure you want to download the \"", 
      filename, "\" fantom5 file?\nFINAL SIZE (might include multiple files); ",
      size, "MB\n" , "Answer \"yes\" or \"y\": "))
    answer <- tolower(readline())
    if (answer == "y"|answer == "yes") {
      print("Starting download!")
    } else {
      print("Download cancelled!")
      # return null to cancel
      return(NULL)
    }
  }
  bfc <- BiocFileCache(ask = FALSE)
  full_filename <- bfcrpath(bfc, url)
  
  if (unpack) {
    # check if unpacked is already in the bioc file cache
    rname <- basename(tools::file_path_sans_ext(full_filename))
    if (!nrow(bfcquery(bfc, query = rname, field = "rname"))) {
    # unpack to bfcnew(bfc, rname) if not available
      gunzip(full_filename, overwrite = TRUE, remove = FALSE)
      bfcnew(bfc, rname)
    }
    full_filename <- tools::file_path_sans_ext(full_filename)
  }
  return(as.character(full_filename))
}
