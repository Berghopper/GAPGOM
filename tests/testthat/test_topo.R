# test lncrnappred functions (all dependent functions)
context("topoicsim prediction functions")

testmaintopo <- function() {
  list1 <- c("126133","221","218","216","8854","220","219","160428","224",
             "222","8659","501","64577","223","217","4329","10840","7915")
  return(
    GAPGOM::topo_ic_sim(list1, list1, ont="MF", organism="human", drop=NULL)
    )
}

test_that("topo algorithms are working", {
  expect_equal(GAPGOM::topo_ic_sim_g1g2("218", "501", ont="MF", 
                                        organism="human", drop=NULL),
               genetopo_baseresult)
  expect_equal(testmaintopo(), maintopo_baseresult)
})