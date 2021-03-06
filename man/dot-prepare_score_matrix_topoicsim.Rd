% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_preparation.R
\name{.prepare_score_matrix_topoicsim}
\alias{.prepare_score_matrix_topoicsim}
\title{GAPGOM internal - .prepare_score_matrix_topoicsim()}
\usage{
.prepare_score_matrix_topoicsim(vec1, vec2, sparse = FALSE,
  old_scores = FALSE)
}
\arguments{
\item{vec1}{vector1 with arbitrary lookup names}

\item{vec2}{vector2 with arbitrary lookup names.}

\item{sparse}{whether to implement Matrix or matrix (default=F --> matrix)}

\item{old_scores}{old score matrix that has overlapping values with 
currently generated score matrix.}
}
\value{
The score matrix with names and NA's.
}
\description{
This function is an internal function and should not be called by the user.
}
\details{
Prepared score dataframe for the semantic measures
}
\section{Notes}{

Internal function used in multiple functions of topoICSim.
}

\keyword{internal}
