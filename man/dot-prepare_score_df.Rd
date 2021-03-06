% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_preparation.R
\name{.prepare_score_df}
\alias{.prepare_score_df}
\title{GAPGOM internal - .prepare_score_df()}
\usage{
.prepare_score_df(original_ids, score, gene_id)
}
\arguments{
\item{original_ids}{rowname ID's for selected expression data.}

\item{score}{score vector}

\item{gene_id}{Gene ID originally selected in the prediction function. 
(gets filtered)}
}
\value{
The score dataframe with ensmbl ID's
}
\description{
This function is an internal function and should not be called by the user.
}
\details{
Prepared score dataframe for the enrichment analysis.
}
\section{Notes}{

Internal function used in prediction_function().
}

\keyword{internal}
