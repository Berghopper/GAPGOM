# test miscellaneous things
# $sessioninfo$otherPkgs$org.Hs.eg.db

context("Miscellaneous tests")

freq_gos_ismaintained <- TRUE

if (freq_gos_ismaintained) {
  
  test_that("freq_go_pairs human db up-to-date", {
    library(org.Hs.eg.db)
    expect_equal(GAPGOM:::freq_go_pairs$
                   sessioninfo$otherPkgs$org.Hs.eg.db$Version,
                 .get_package_version("org.Hs.eg.db"))
  })
  
  test_that("freq_go_pairs mouse db up-to-date", {
    library(org.Mm.eg.db)
    expect_equal(GAPGOM:::freq_go_pairs$
                   sessioninfo$otherPkgs$org.Mm.eg.db$Version, 
                 .get_package_version("org.Mm.eg.db"))
  })
}
