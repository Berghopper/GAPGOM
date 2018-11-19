# test lncrnappred functions (all dependent functions)
context("topoicsim prediction functions")

testmaintopo <- function() {
  list1 <- c("126133","221","218","216","8854","220","219","160428","224",
             "222","8659","501","64577","223","217","4329","10840","7915")
  return(
    GAPGOM::topo_ic_sim_genes("human", "MF", list1, list1, drop=NULL)
    )
}

test_that("topo algorithms are working", {
  expect_equal(GAPGOM::topo_ic_sim_term("human", "MF", "GO:0018478", "GO:0047105"), 
               termtopo_baseresult)
  expect_equal(GAPGOM::topo_ic_sim_genes("human", "MF", "218", "501", 
                                         drop=NULL),
               genetopo_baseresult)
  expect_equal(testmaintopo(), maintopo_baseresult)
})