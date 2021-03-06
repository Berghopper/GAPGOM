#' GAPGOM - Expression Data.
#'
#' A Bioconductor ExpressionSet object containing a subset from the 
#' lncRNA2Function data.[1]
#' Because the original lncRNA2Function data was unavailable for some time,
#' this subset is selected with expression data from (this is based on the
#' lncRNA2Function however): https://tare.medisin.ntnu.no/pred_lncRNA/ .
#' The original lncRNA2Function data is now online on a new domain:
#' http://bio-annotation.cn/lncrna2function/
#'
#' @format An ExpressionSet containing 1001 rows of the lncRNA2Function 
#' expression dataset. 
#' Annotation data can be found under;
#' pData(featureData(GAPGOM::expset))
#' and expression values under;
#' assayData(GAPGOM::expset)[["exprs"]]
#'
#' @source \url{http://bio-annotation.cn/lncrna2function/}
#' @references [1]. Jiang Q et al.: \strong{LncRNA2Function: a comprehensive 
#' resource for functional investigation of human lncRNAs based on RNA-seq 
#' data.} in: \emph{BMC Genomics}, 2015, doi: 10.1186/1471-2164-16-S3-S2.
"expset"