## ----setup, include = FALSE----------------------------------------------
library(knitr)
library(kableExtra)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  error = F
)

## ----f5, eval=F----------------------------------------------------------
#  # download the fantom5 data file
#  fantom_file <- fantom_download("./", organism = "mouse",
#                                 noprompt = T) # saves filename
#  # load the file (use fantom_file variable if doing all at once)
#  ft5 <- fantom_load_raw("./mm9.cage_peak_phase1and2combined_tpm_ann.osc.txt",
#  verbose = T)
#  # remove first two rows from fantom5 data (these are seperate statistis,
#  # we just need expressionvalues)
#  ft5$df <- ft5$df[3:nrow(ft5$df),]
#  
#  # convert the raw fantom table to an ExpressionSet
#  expset <- fantom_to_expset(ft5, verbose = T)

## ----randvals, eval=F----------------------------------------------------
#  # select 10000 random IDs
#  go_data <- GAPGOM::set_go_data("human", "BP", computeIC = F)
#  random_ids <- unique(sample(go_data@geneAnno$ENTREZID, 10000)) # and only keep
#  # uniques
#  
#  # make general dataframe.
#  expressions <- data.frame(random_ids)
#  colnames(expressions) <- "ENTREZID"
#  expressions$ID
#  
#  # n expression values depending on the amount of unique IDs that are present
#  expressionvalues <- abs(rnorm(length(random_ids)*6))*10000
#  expressions[,2:7] <- expressionvalues
#  View(expressions)

## ----expset, eval=F------------------------------------------------------
#  expression_matrix <- as.matrix(expressions[,2:ncol(expressions)])
#  rownames(expression_matrix) <- expressions$ENTREZID
#  featuredat <- as.data.frame(expressions$ENTREZID) # And everything else besides expressionvalues (preferably you don't even need to include the IDs themselves here!)
#  rownames(featuredat) <- expressions$ENTREZID # because they will be the rownames anyway.
#  expset <- ExpressionSet(expression_matrix,
#                          featureData = new("AnnotatedDataFrame",
#                          data=featuredat))
#  
#  # To see how it is structured;
#  View(expset)
#  View(expset@assayData$exprs) # where expressionvalues are stored.
#  View(expset@featureData@data) # where other information is stored.

