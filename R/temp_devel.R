# DEVELOPMENT R SCRIPTING FILE - TEMPORARY
library(biomaRt)
library(CAGEr)


#' ID convert script
#' @import biomaRt
ens_to_entrez <- function(id) {
  ensembl = useMart("ensembl",dataset="hsapiens_gene_ensembl")
}

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
  # '/media/casper/USB_ccpeters/internship_thesis/data/hg19.cage_peak_phase1and2combined_counts.osc.txt'
  print("reading columnvariables...")
  con <- file(filepath, "r")
  linecount <- 0
  while ( TRUE ) {
    line <- readLines(con, n = 1)
    if (!startsWith(line, "#")) {
      break
    }
    linecount <- linecount + 1
  }
  close(con)
  return(fread(filepath, skip=linecount))
}

fantom_load <- function() {
  astrocyteSamples <- FANTOM5humanSamples[grep("Astrocyte", 
                                               FANTOM5humanSamples[,"description"]),]
  return (importPublicData(source = "FANTOM5", dataset = "human", 
                           sample = astrocyteSamples[1:3,"sample"]))
} 

# try correct to RPKM vals
ft5 <- fantom_load()
View(ft5)
bla <- apply(ft5@tagCountMatrix, 1, function(x) {x/(nrow(ft5@tagCountMatrix)/1000000)})
