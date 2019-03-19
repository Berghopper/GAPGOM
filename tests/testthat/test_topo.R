# test lncrnappred functions (all dependent functions)
context("topoicsim prediction functions")

test_that("topoicsim term level", {
  expect_equal(GAPGOM::topo_ic_sim_term("human", "MF", "GO:0018478", 
                                        "GO:0047105"), 
               gapgom_tests$termtopo_baseresult)
})

test_that("topoicsim gene level", {
  #skip("R cmd check time constraints")
  expect_equal(GAPGOM::topo_ic_sim_genes("human", "MF", "218", "501",
                                         use_precalculation = FALSE,
                                         drop=NULL, progress_bar = FALSE),
               gapgom_tests$genetopo_baseresult)
})

testgenesettopo <- function() {
  list1 <- c("126133","221","218","216","8854","220","219","160428","224",
             "222","8659","501","64577","223","217","4329","10840","7915",
             "5832")
  return(
    GAPGOM::topo_ic_sim_genes("human", "MF", list1[1:5], list1[1:5],
                              use_precalculation = TRUE, 
                              progress_bar = FALSE)
  )
}

test_that("topoicsim geneset level", {
  #skip("R cmd check time constraints")
  expect_equal(testgenesettopo(), gapgom_tests$maintopo_baseresult)
})
