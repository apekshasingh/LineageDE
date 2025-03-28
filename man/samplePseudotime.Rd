% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/samplePseudotime.R
\name{samplePseudotime}
\alias{samplePseudotime}
\title{samplePseudotime}
\usage{
samplePseudotime(
  cellData,
  conditions,
  clusters,
  start.cluster,
  end.clusters,
  experiments = NULL,
  nSamples = 100,
  pseudo.cond = NULL,
  pseudo.data = NULL,
  parallel = FALSE,
  n_cores = NULL,
  seed = 123
)
}
\arguments{
\item{cellData}{matrix (size = cells by dimensions) of PCA projection (or other dim reduc) representation of cells}

\item{conditions}{factor (length = cells) specifying condition to which each cell belongs}

\item{clusters}{factors (length = cells) specifying to which input cluster each cell belongs.  Input clusters allow slingshot to create an initial MST.}

\item{start.cluster}{string specifying the name of the cluster from which lineages should begin}

\item{end.clusters}{vectors of string specifying names of clusters where lineages should terminate}

\item{experiments}{factor (length = cells) specifying replicate/batch to which each cell belongs if applicable (otherwise NULL)}

\item{nSamples}{integer specifying number of samples to generate (defaults to 100)}

\item{pseudo.cond}{specifies name of condition to be used to build pseudotime trajectory (defaults to NULL and uses all cells)}

\item{pseudo.data}{matrix of alternative cells' PCA projections (or other dim reduc) to use to build pseudotime trajectories.  cellData is then projected onto these trajectories.}

\item{parallel}{boolean specifying whether or not to run in parallel (requires BiocParallel)}

\item{n_cores}{integer specifying number of cores to use if running in parallel}

\item{seed}{set random seed}
}
\value{
list of pseudotime inputs utilized for LineageDE (ori.pseudotimes, ori.weights, samples, sub.pseudotimes, sub.weights)
}
\description{
Provides functionality to generate multiple samplings of pseudotime using the package
slingshot to generate the pseudotime trajectories.  This output can be used for input
with LineageDE.  This function can be modified to incorporate other pseudotime algorithms.
}
