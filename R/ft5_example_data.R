#' GAPGOM - FANTOM5 Expression data.
#'
#' A dataframe containing fantom5 expression values from the mouse genus.
#' It does not contain all values as this would be too much to include with
#' the package. [1]
#'
#' The host url/website is unavailable.
#'
#' @format An ExpressionSet containing the first 15 columns and 10000 rows of
#' the mouse RLE normalized fantom5 expression dataset where entrezids are
#' present. Annotation data can be found under;
#' GAPGOM::ft5_example_data@featureData@data
#' and expression values under;
#' GAPGOM::ft5_example_data@assayData$exprs
#'
#' @source \url{http://fantom.gsc.riken.jp/5/data/}
#' @references [1]. Derek de Rie, Imad Abugessaisa et al.: \strong{An 
#' integrated expression atlas of miRNAs and their promoters in human and mouse}
#' \emph{Nature 2017}.
"ft5_example_data"