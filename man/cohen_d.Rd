% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Desc_function.R
\name{cohen_d}
\alias{cohen_d}
\title{cohen_d}
\usage{
cohen_d(Data, index_1, index_2)
}
\arguments{
\item{Data}{A data frame that contains sample sizes, means, and standard
deviations.}

\item{index_1}{Number of group 1.}

\item{index_2}{Number of group_2}
}
\value{
A vector of pairwise Cohen'ds.
}
\description{
This function estimates Cohen's d in \code{\link{MAGMA_desc}}.
}
\details{
Inner function of \code{\link{MAGMA_desc}} that computes Cohen's d using the
pooled SD.
}
\author{
Julian Urban
}
