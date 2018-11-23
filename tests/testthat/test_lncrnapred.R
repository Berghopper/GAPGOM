# test lncrnappred functions (all dependent functions)
context("lncRNA prediction functions")

lncrnapred_test <- function() {
   # Example with default dataset, take a look at the data documentation
   # to fully grasp what's going on with making of the filter etc. (Biobase 
   # ExpressionSet)
   
   # keep everything that is a protein coding gene
   filter_vector <- expset@featureData@data[(expset@featureData@data$GeneType=="protein_coding"),]$GeneID
   # set gid and run.
   gid <- "ENSG00000228630"
   
   return(GAPGOM::expression_prediction(gid, 
                                        GAPGOM::expset, 
                                        "human", 
                                        "BP",
                                        id_translation_df = GAPGOM::id_translation_df,
                                        id_select_vector = filter_vector,
                                        method = "combine", 
                                        verbose = TRUE, 
                                        filter_pvals = TRUE)
   )
}

test_that("expression_prediction_function is working", {
  expect_equal(lncrnapred_test(), lncrnapred_baseresult)
})
