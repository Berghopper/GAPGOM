## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval=F-------------------------------------------------------------
#  # select 10000 random IDs
#  go_data <- GAPGOM:::.set_go_data("human", "BP", computeIC = F)
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

## ---- eval=F-------------------------------------------------------------
#  expression_matrix <- as.matrix(expressions[,2:ncol(expressions)])
#  rownames(expression_matrix) <- expressions$ENTREZID
#  featuredat <- as.data.frame(expressions$ENTREZID) # And everything else besides
#  # expressionvalues (preferably you don't even need to include the IDs themselves
#  # here!)
#  rownames(featuredat) <- expressions$ENTREZID # because they will be the rownames
#  # anyway.
#  expset <- ExpressionSet(expression_matrix,
#                          featureData = new("AnnotatedDataFrame",
#                          data=featuredat))
#  
#  # To see how it is structured;
#  View(expset)
#  View(expset@assayData$exprs) # where expressionvalues are stored.
#  View(expset@featureData@data) # where other information is stored.

## ------------------------------------------------------------------------
# Example with default dataset, take a look at the data documentation
# to fully grasp what's going on with making of the filter etc. (Biobase 
# ExpressionSet)

# make a filter of all the ids that contain a 0 in their expression row.
sort_list <- rownames(GAPGOM::ft5_example_data@assayData$exprs[
apply(GAPGOM::ft5_example_data@assayData$exprs, 1, 
function(row) {all(row==0)}),])
# set an artbitrary gene you want to find similarities for. (name of 5th row 
# in this case)
gid <- rownames(GAPGOM::ft5_example_data@assayData$exprs)[5]
# get the result.
result <- GAPGOM::expression_prediction(gid, 
                                       GAPGOM::ft5_example_data, 
                                       sort_list, 
                                       "mouse", 
                                       "BP"
                                       )
result

## ------------------------------------------------------------------------
# single gene mode
result <- GAPGOM::topo_ic_sim_genes("218", "501", ont="MF", organism="human", 
                                   drop=NULL)

result
# genelist mode
list1 <- c("126133","221","218","216","8854","220","219","160428","224",
"222","8659","501","64577","223","217","4329","10840","7915")
result <- GAPGOM::topo_ic_sim_genes(list1, list1, ont="MF", organism="human", 
                              drop=NULL)
result

# calculate intersetsim (intra if gene lists would be different)
mean(result$GeneSim)

## ------------------------------------------------------------------------
sessionInfo()

