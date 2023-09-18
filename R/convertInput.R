#' convertInput
#'
#' To compare gene expression between two or more lineages from the same pseudotime analysis
#' rather than comparing gene expression between conditions across single lineage this function
#' can be utilized to convert the pseudotime input to make it compatible with LineageDE.
#' A copy of each cell is created for each lineage and the output conditions will
#' correspond to the different lineages.
#'
#' @param counts matrix (size = cells by genes) of raw (integer) count values of gene expression
#' @param ori.pseudotimes matrix (size = cells by lineages) of pseudotime values
#' @param linsCompare vector of lineage indices to be compared (e.g. to compare lineage 1 to lineage 2, linsCompare = c(1, 2))
#' @param ori.weights matrix (size = cells by lineages) of weights from pseudotime analysis if applicable (otherwise NULL)
#' @param experiments factor (length = cells) specifying replicate/batch to which each cell belongs if applicable (otherwise NULL)
#' @param samples matrix (size = cells by num_samples) specifying indices cells utilized in each sample of pseudotime if applicable (otherwise NULL)
#' @param sub.pseudotimes matrix (size = cells by lineages by num_samples) samples of pseudotime values if applicable (otherwise NULL)
#' @param sub.weights matrix (size = cells by lineages by num_samples) samples of pseudotime weights if applicable (otherwise NULL)
#' @param seed set random seed
#'
#' @returns converted input in the form of a list containing counts, pseudotimes, weights, experiments (if applicable), samples (if applicable), and the generated conditions factor to use with LineageDE
#'
#' @export
convertInput<- function(counts, ori.pseudotimes, linsCompare, ori.weights=NULL, experiments=NULL, samples=NULL, sub.pseudotimes=NULL, sub.weights=NULL, seed=123){

  #To compare lineages within the same experiment we convert each lineage to a unique condition by replicating cells

  ## DATA INPUTS:
  # assuming multiple subsamplings of pseudotime ouptut already generated

  # counts: matrix (cell by genes - raw values)

  # pseudotime results:
  # ori.pseudotimes = cells x lineages
  # ori.weights = cells x lineages
  # samples = indices (80% cells) x n_samples
  # sub.pseudotimes = cells (80%) x lineages x n_samples
  # sub.weights = cells (80%) x lineages x n_samples
  # if no weights provided, assume all cells have equal weight
  if (is.null(ori.weights)){ori.weights = array(1, dim(ori.pseudotimes))}
  if (is.null(sub.weights) & !is.null(samples)){sub.weights = array(1, dim(sub.pseudotimes))}
  # all pseudotime samples should be normalized! (max 1)

  # linsCompare provides the indices of the lineages to be compared
  ori.pseudotimes = ori.pseudotimes[, linsCompare]
  ori.weights = ori.weights[, linsCompare]
  if (!is.null(samples)){
  sub.pseudotimes = sub.pseudotimes[, linsCompare, ]
  sub.weights = sub.weights[, linsCompare, ]}

  # Convert input for LineageDE

  nLin = ncol(ori.pseudotimes)
  nCell = dim(counts)[1]
  lin_assign = rep(seq_len(nLin), nCell)
  lin_assign = factor(lin_assign)
  counts = counts[rep(1:nCell, each=nLin), ]
  if (!is.null(experiments)){
  experiments = experiments[rep(1:nCell, each=nLin)]}

  collapseMatrix<-function(mat, nCell, nLin){
    mat<-t(mat)
    dim(mat)<-c(nCell*nLin, 1)
    mat<-as.matrix(mat, ncol=1)
    return(mat)
  }
  ori.pseudotimes = collapseMatrix(ori.pseudotimes, nCell, nLin)
  ori.weights = collapseMatrix(ori.weights, nCell, nLin)

  if (!is.null(samples)){
  sample.nCell = dim(samples)[1]
  nSamples = ncol(samples)

  convertIndex<-function(sample_indices){unlist(lapply(sample_indices, function(index){((nLin*(index-1))+1):((nLin*(index-1))+nLin)}))}
  samples<-apply(samples, 2, convertIndex)

  sub.pseudotimes_temp = array(dim = c(sample.nCell*nLin, 1, nSamples))
  sub.weights_temp = array(dim = c(sample.nCell*nLin, 1, nSamples))


  for (s in 1:nSamples){
    sub.pseudotimes_temp[, , s]<-collapseMatrix(sub.pseudotimes[, , s], sample.nCell, nLin)
    sub.weights_temp[, , s]<-collapseMatrix(sub.weights[, , s],sample.nCell, nLin)
  }

  sub.pseudotimes = sub.pseudotimes_temp
  rm(sub.pseudotimes_temp)

  sub.weights = sub.weights_temp
  rm(sub.weights_temp)}

  newInput<-list(counts=counts, ori.pseudotimes=ori.pseudotimes, ori.weights=ori.weights, conditions=lin_assign, experiments=experiments, samples=samples, sub.pseudotimes=sub.pseudotimes, sub.weights=sub.weights)
  return(newInput)
}
