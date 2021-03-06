% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data_preparation.R
\name{.set_identical_items}
\alias{.set_identical_items}
\title{GAPGOM internal - .set_identical_items()}
\usage{
.set_identical_items(score_matrix)
}
\arguments{
\item{score_matrix}{score matrix for topoclsim}
}
\value{
Same matrix with correct sets set to 1.
}
\description{
This function is an internal function and should not be called by the user.
}
\details{
Sets identical items within a score matrix, the similarity will always be 1.
}
\section{Notes}{

Internal function used in various topoclsim algorithms.
}

\keyword{internal}
