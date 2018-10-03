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

# strsplit(ft5_ranncom$df[!(is.na(ft5_ranncom$df$`entrezgene (genes) id associated with the transcript`) | ft5_ranncom$df$`entrezgene (genes) id associated with the transcript`==""),][1,5][[1]], "\\:")[[1]][2]

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

#' convert ft5 to expset
#' @examples 
#' ft5_r <- fantom_load_raw("/media/casper/USB_ccpeters/internship_thesis/data/f5/mouse/mm9.cage_peak_phase1and2combined_tpm_ann.osc.txt", verbose=T)
#' fantom_to_expset(ft5_r)
#' @importFrom Biobase ExpressionSet annotatedDataFrameFrom
fantom_to_expset <- function(fanraw, verbose = F) {
  fan <- fanraw$df
  meta <- fanraw$meta
  if (verbose) print("filtering out empty entrez ids...")
  fan <- fan[!(is.na(fan$`entrezgene (genes) id associated with the transcript`) | fan$`entrezgene (genes) id associated with the transcript`==""),]
  #fan <- fan[1:10000,1:15]
  if (verbose) print("DONE")
  if (verbose) print("converting to expressionset...")
  expression_matrix <- as.matrix(fan[,7:ncol(fan)])
  rownames(expression_matrix) <- fan$`CAGE peak id`
  featuredat <- as.data.frame(fan[,2:6])
  rownames(featuredat) <- fan$`CAGE peak id`
  print(nrow(featuredat))
  print(nrow(expression_matrix))
  return(ExpressionSet(expression_matrix, featureData = new("AnnotatedDataFrame", data=featuredat),...=meta))
  
}

#' @importFrom GEOquery gunzip
download_fantom5 <- function(down_dir, organism="human") {
  url <- switch(organism, "human" = "http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz",
         "mouse" = "http://fantom.gsc.riken.jp/5/datafiles/latest/extra/CAGE_peaks/mm9.cage_peak_phase1and2combined_tpm_ann.osc.txt.gz")
  filename <- strsplit(url, "/")[[1]][length(strsplit(url, "/")[[1]])]
  full_filename <- paste0(down_dir,"/",filename)
  download.file(url, full_filename, "auto")
  gunzip(full_filename, overwrite = T)
}



#' convert into expressionset
fantom_load_expressionset <- function(filepath, colselector, verbose = F) {
  fanraw <- fantom_load_raw(filepath, verbose = verbose)
  return(fantom_to_expset(fanraw, colselector, verbose = verbose))
}

#' @export
fantom_test_data <- function() {
  ft5 <- fantom_load_raw("/media/casper/USB_ccpeters/internship_thesis/data/f5/mouse/mm9.cage_peak_phase1and2combined_tpm_ann.osc.txt", verbose=T)
  return(fantom_to_expset(ft5))
}


