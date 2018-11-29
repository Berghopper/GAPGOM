# test functions from parsing files

context("parsing functions")

test_that("entrez raw (fantom) id to entrez", {
  expect_equal(.entrezraw_to_entrez("entrezgene:123"), "123")
  expect_equal(.entrezraw_to_entrez(c("entrezgene:123", "entrezgene:345")), c("123", "345"))
})

test_that("GO quantification", {
  testdf <- data.frame(list(GO=c("GO:123", "GO:123", "GO:1234")))
  testdf2 <- data.frame(list(FAKEID=c("1", "2", "3"), 
                             GO=c("GO:123", "GO:123", "GO:1234")))
  result_structure <- structure(list(
    GO = structure(1:2, .Label = c("GO:123", "GO:1234"
  ), class = "factor"), N = 2:1), row.names = c(NA, -2L), class = "data.frame")
  expect_equal(.term_id_to_ext_id(testdf), result_structure)
  expect_equal(.term_id_to_ext_id(testdf2), result_structure)
})

test_that("GO quantification subset", {
  expect_equal(.ext_id_to_term_id(id_translation_df[
    id_translation_df$ORIGID %in% unique(id_translation_df$ORIGID)[1:10],], 
    c("ENSG00000148516", "ENSG00000006534")), gapgom_tests$goquant2)
})
