# test lncrnappred functions (all dependent functions)
context("lncRNA prediction functions")

lncrnapred_test <- function() {
  # Example with default dataset, take a look at the data documentation
  # to fully grasp what's going on with making of the filter etc. (Biobase 
  # ExpressionSet)
  
  # make a filter of all the ids that contain a 0 in their expression row.
  sort_list <- rownames(GAPGOM::ft5_example_data@assayData$exprs[
    apply(GAPGOM::ft5_example_data@assayData$exprs, 1, 
          function(row) {all(row==0)}),])
  # set an artbitrary gene you want to find similarities for. (5th row in this
  # case)
  gid <- rownames(ft5_example_data@assayData$exprs)[5]
  return(GAPGOM::expression_prediction(gid, 
                                         GAPGOM::ft5_example_data, 
                                         sort_list, 
                                         "mouse", 
                                         "BP"
  ))
}

test_that("expression_prediction_function is working", {
  options(stringsAsFactors = F)
  expect_equal(lncrnapred_test(), lncrnapred_baseresult)
})