#' GAPGOM - Expression Data.
#'
#' An ExpressionSet from Bioconductor containing a subset from the 
#' lncRNA2Function data.[1]
#' The original lncRNA2Function data sadly seems permanently unavailable,
#' however, an intermediary selection of the dataset is still available.
#'
#' https://tare.medisin.ntnu.no/pred_lncRNA/
#'
#' @format An ExpressionSet containing 1001 rows of the lncRNA2Function 
#' expression dataset. 
#' Annotation data can be found under;
#' pData(featureData(GAPGOM::expset))
#' and expression values under;
#' assayData(GAPGOM::expset)[["exprs"]]
#'
#' @source \url{https://tare.medisin.ntnu.no/pred_lncRNA/}
#' @references [1]. Jiang Q et al.: \strong{LncRNA2Function: a comprehensive 
#' resource for functional investigation of human lncRNAs based on RNA-seq 
#' data.} in: \emph{BMC Genomics}, 2015, doi: 10.1186/1471-2164-16-S3-S2.
"expset"