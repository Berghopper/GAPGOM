#' GAPGOM -  EnsemblID -> GOID conversion table
#'
#' A dataframe containing EnsemblIDs and corresponding GOIDs for refering back
#' and forth between the IDs.
#'
#' Source unknown.
#'
#' @format A data frame with 252887 rows and 3 variables/columns:
#' \describe{
#'   \item{ensembl_gene_id}{The EnsemblID of the gene.}
#'   \item{go_id}{The GOID of the property. Multiple GOIDs per EnsemblIDs
#'   are possible.}
#'   \item{namespace_1003}{Toplevel GO node.}
#' }
#'
"ensembl_id_to_go_id"
