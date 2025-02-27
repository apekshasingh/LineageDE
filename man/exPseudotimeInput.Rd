% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{exPseudotimeInput}
\alias{exPseudotimeInput}
\title{Example Input for samplePsuedotime}
\usage{
data(exPseudotimeInput)
}
\description{
A list containing example input data to run samplePseudotime_HSPC
}
\details{
\describe{
  \item{cellData}{matrix of PCA projection (cells x dimensions)}
  \item{conditions}{factor specifying experimental conditions}
  \item{experiments}{factor specifying experimental replicates}
  \item{clusters}{factor specifying input cluster for each cell (input for Slingshot)}
}
}
\keyword{datasets}
