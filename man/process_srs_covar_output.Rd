% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_srs_covar_output.R
\name{process_srs_covar_output}
\alias{process_srs_covar_output}
\title{Process the output of srs_covar for PCA, etc.}
\usage{
process_srs_covar_output(A, npp = 5)
}
\arguments{
\item{A}{the output list from \code{srs_covar()}}

\item{npp}{the number of principal components to retain in a tibble for plotting
and facetting over the different PCs.}
}
\description{
Compute a lot of stuff and returns a list
}
