# DEVELOPMENT R SCRIPTING FILE - TEMPORARY

library(CAGEr)
library(biomaRt)

# temporary conversion for old to new data.

# expr matrix
new_exp <- as.matrix(GAPGOM::expression_data[,4:22])
rownames(new_exp) <- GAPGOM::expression_data[,1]
# feature df
fdat <- GAPGOM::expression_data[,2:3]
rownames(fdat) <- GAPGOM::expression_data[,1]
fdat_andf <- as(fdat, "AnnotatedDataFrame")

# new ExpressionSet object.
new_obj <- Biobase::ExpressionSet(new_exp, featureData = fdat_andf)

# fantom5 data load.
fantom_load_raw <- function(filepath) {
  options(stringsAsFactors = F)
  library(data.table)
  # '/media/casper/USB_ccpeters/internship_thesis/data/hg19.cage_peak_phase1and2combined_counts.osc.txt'
  print("reading columnvariables...")
  con <- file(filepath, "r")
  translation_df <- data.frame(c(0,0,0))
  linecount <- 0
  while ( TRUE ) {
    line <- readLines(con, n = 1)
    if (!startsWith(line, "#")) {
      break
    } else {
      filtered_line <- sapply(strsplit(line, "##|\\[|\\]|\\="), function(x){x[!x ==""]})
      translation_df <- data.frame(translation_df, filtered_line)
    }
    linecount <- linecount + 1
  }
  translation_df <- as.data.frame(t(translation_df))
  rownames(translation_df) <- 1:nrow(translation_df)
  colnames(translation_df) <- c("type","key","value")
  translation_df <<- translation_df
  close(con)
  fan <<- fread(filepath, skip=linecount)
  colnames(fan) <- sapply(colnames(fan), function(x){translation_df[translation_df$key==x,]$value})
  # now make a ExpressionSet LEFT OFF HERE
  translation_df[translation_df$type!="ColumnVariables" | translation_df$type!=0,]
  return(fan)
}

ft5_r <- fantom_load_raw("/media/casper/USB_ccpeters/internship_thesis/data/hg19.cage_peak_phase1and2combined_counts.osc.txt")
ft5_rann <- fantom_load_raw("/media/casper/USB_ccpeters/internship_thesis/data/hg19.cage_peak_phase1and2combined_tpm_ann.osc.txt")

fantom_load <- function() {
  astrocyteSamples <- FANTOM5humanSamples[grep("Astrocyte", 
                                               FANTOM5humanSamples[,"description"]),]
  return (importPublicData(source = "FANTOM5", dataset = "human", 
                           sample = astrocyteSamples[1:6,"sample"]))
} 

# try correct to RPKM vals
ft5 <- fantom_load()
View(ft5)
bla <- t(apply(ft5@tagCountMatrix, 1, function(x) {x/(nrow(ft5@tagCountMatrix)/1000000)}))
rownames(bla) <- with(ft5@CTSScoordinates, paste(chr, pos, strand, sep="_"))
bla
