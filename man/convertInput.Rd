% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/convertInput.R
\name{convertInput}
\alias{convertInput}
\title{convertInput}
\usage{
convertInput(
  counts,
  ori.pseudotimes,
  linsCompare,
  ori.weights = NULL,
  experiments = NULL,
  samples = NULL,
  sub.pseudotimes = NULL,
  sub.weights = NULL,
  seed = 123
)
}
\arguments{
\item{counts}{matrix (size = cells by genes) of raw (integer) count values of gene expression}

\item{ori.pseudotimes}{matrix (size = cells by lineages) of pseudotime values}

\item{linsCompare}{vector of lineage indices to be compared (e.g. to compare lineage 1 to lineage 2, linsCompare = c(1, 2))}

\item{ori.weights}{matrix (size = cells by lineages) of weights from pseudotime analysis if applicable (otherwise NULL)}

\item{experiments}{factor (length = cells) specifying replicate/batch to which each cell belongs if applicable (otherwise NULL)}

\item{samples}{matrix (size = cells by num_samples) specifying indices cells utilized in each sample of pseudotime if applicable (otherwise NULL)}

\item{sub.pseudotimes}{matrix (size = cells by lineages by num_samples) samples of pseudotime values if applicable (otherwise NULL)}

\item{sub.weights}{matrix (size = cells by lineages by num_samples) samples of pseudotime weights if applicable (otherwise NULL)}

\item{seed}{set random seed}
}
\value{
converted input in the form of a list containing counts, pseudotimes, weights, experiments (if applicable), samples (if applicable), and the generated conditions factor to use with LineageDE
}
\description{
To compare gene expression between two or more lineages from the same pseudotime analysis
rather than comparing gene expression between conditions across single lineage this function
can be utilized to convert the pseudotime input to make it compatible with LineageDE.
A copy of each cell is created for each lineage and the output conditions will
correspond to the different lineages.
}
