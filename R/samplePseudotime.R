#' samplePseudotime
#'
#' Provides functionality to generate multiple samplings of pseudotime using the package
#' slingshot to generate the pseudotime trajectories.  This output can be used for input
#' with LineageDE.  This function can be modified to incorporate other pseudotime algorithms.
#'
#' @param cellData matrix (size = cells by dimensions) of PCA projection (or other dim reduc) representation of cells
#' @param conditions factor (length = cells) specifying condition to which each cell belongs
#' @param clusters factors (length = cells) specifying to which input cluster each cell belongs.  Input clusters allow slingshot to create an initial MST.
#' @param start.cluster string specifying the name of the cluster from which lineages should begin
#' @param end.clusters vectors of string specifying names of clusters where lineages should terminate
#' @param experiments factor (length = cells) specifying replicate/batch to which each cell belongs if applicable (otherwise NULL)
#' @param nSamples integer specifying number of samples to generate (defaults to 100)
#' @param pseudo.cond specifies name of condition to be used to build pseudotime trajectory (defaults to NULL and uses all cells)
#' @param pseudo.data matrix of alternative cells' PCA projections (or other dim reduc) to use to build pseudotime trajectories.  cellData is then projected onto these trajectories.
#' @param parallel boolean specifying whether or not to run in parallel (requires BiocParallel)
#' @param n_cores integer specifying number of cores to use if running in parallel
#' @param seed set random seed
#'
#' @returns list of pseudotime inputs utilized for LineageDE (ori.pseudotimes, ori.weights, samples, sub.pseudotimes, sub.weights)
#'
#' @export
samplePseudotime<- function(cellData, conditions, clusters, start.cluster, end.clusters, experiments=NULL, nSamples=100, pseudo.cond=NULL, pseudo.data = NULL, parallel=FALSE, n_cores=NULL, seed=123){

  #cellData is the PCA projection (or other dimensionality reduced representation) of the gene expression dataset
  #conditions is a factor specifying the condition to which each cell belongs
  #experiments is a factor specifying the experiment/batch to which each cell belongs

  #clusters is a factor specifying the cluster to which each cell belongs (corresponds to input for Slingshot)
  #start.cluster, end.clusters specify which clusters mark the beginning and ends of lineages
  #nSamples is the number of pseudotime samples to generate

  #pseudo.cond specifies the condition from which the pseudotime trajectories will be built,
  # the other conditions are projected onto these trajectories
  # (if all conditions should be utilized, simply leave as null)

  #pseudo.data is the PCA projection (or other dimensionality reduced representation)
  # of an alternative dataset from which the pseudotime trajectories are constructed
  # and cellData is projected onto (leave as null if cellData should be used)

  #seed allows random samples to be consistent across multiple function calls
  set.seed(seed)

  clusters <-droplevels(clusters)
  conditions <-droplevels(conditions)
  if(is.null(experiments)){experiments = factor(rep("expt", length(conditions)))
  } else {experiments <- droplevels(experiments)}

  getSlingshotResults<-function(cellData, conditions, clusters, start.cluster, end.clusters, pseudo.cond, pseudo.data){

    nCell= nrow(cellData)
    nLin = length(end.clusters)
    pseudotimes = array(dim = c(nCell, nLin))
    weights = array(dim = c(nCell, nLin))
    distances = array(dim = c(nCell, nLin))

    if (is.null(pseudo.cond) & is.null(pseudo.data)){
      slingshot_results<-slingshot::slingshot(cellData, clusters, start.clus = start.cluster, end.clus=end.clusters)
      max_pt = max(slingshot::slingPseudotime(slingshot_results), na.rm=TRUE)
    }
    else if (!is.null(pseudo.cond)){
      ref_cells = conditions==pseudo.cond
      #obtain pseudotime reference
      slingshot_results<-slingshot::slingshot(cellData[ref_cells, ], clusters[ref_cells], start.clus = start.cluster, end.clus=end.clusters)
      max_pt = max(slingshot::slingPseudotime(slingshot_results), na.rm=TRUE)
      #project cells onto reference
      slingshot_results<-slingshot::predict(slingshot_results, cellData)
    }
    else{
      #obtain pseudotime reference (alternative dataset)
      slingshot_results<-slingshot::slingshot(pseudo.data, clusters, start.clus = start.cluster, end.clus=end.clusters)
      max_pt = max(slingshot::slingPseudotime(slingshot_results), na.rm=TRUE)
      #project all cells onto reference
      slingshot_results<-slingshot::predict(slingshot_results, rbind(cellData, pseudo.data))}

    #lineages = slingshot_results@lineages
    #curves = slingshot_results@curves

    lineages = slingshot_results@metadata$lineages
    curves = slingshot_results@metadata$curves

    for (l in 1:nLin){
      lin_id = utils::tail(lineages[[l]], 1)
      lin_id = which(lin_id == end.clusters)

      pseudotimes[, lin_id] = curves[[l]]$lambda[1:nCell]
      weights[, lin_id] = curves[[l]]$w[1:nCell]}

    #max-normalize pseudotimes
    pseudotimes = pseudotimes/max_pt
    return(list(pseudotimes=pseudotimes, weights=weights))
  }

  print("running original pseudotime analysis")
  slingshot_results<-getSlingshotResults(cellData, conditions, clusters, start.cluster, end.clusters, pseudo.cond, pseudo.data)
  ori.pseudotimes <- slingshot_results$pseudotimes
  ori.weights <- slingshot_results$weights

  nCell= nrow(cellData)
  nLin = length(end.clusters)

  print("running sampling pseudotime analysis")

  #if need to sample without replacement
  #samples = matrix(nrow = as.integer(0.8*nCell), ncol=nSamples)
  #sub.pseudotimes = array(dim = c(as.integer(0.8*nCell), nLin, nSamples))
  #sub.weights = array(dim = c(as.integer(0.8*nCell), nLin, nSamples))

  samples = matrix(nrow = nCell, ncol=nSamples)
  sub.pseudotimes = array(dim = c(nCell, nLin, nSamples))
  sub.weights = array(dim = c(nCell, nLin, nSamples))

  getSampling<-function(n, cellData, conditions, clusters, start.cluster, end.clusters, pseudo.cond, pseudo.data){
    print(paste("generating sample", n))

    # Simple Sampling
    #if need to sample without replacement
    #sample_n <- sample(nCell, as.integer(0.8*nCell))
    #else
    #sample_n <- sample(nCell, nCell, replace=TRUE)

    # Stratified Sampling
    sample_n <- c()
    if (is.null(pseudo.data)){
    for (clus in levels(clusters)){
      for (cond in levels(conditions)){
        for (expt in levels(experiments)){
        strat = (clusters==clus) & (conditions == cond) & (experiments==expt)
        #if need to sample without replacement
        #sample_n = append(sample_n, sample(which(strat), as.integer(0.8*sum(strat)))
        sample_n = append(sample_n, sample(which(strat), sum(strat), replace=TRUE))
        }}}
      slingshot_results_n <- getSlingshotResults(cellData[sample_n, ], conditions[sample_n], clusters[sample_n], start.cluster, end.clusters, pseudo.cond, pseudo.data)}
    else {
      sample_p = c()
    for (cond in levels(conditions)){
      for (expt in levels(experiments)){
        strat = (conditions == cond) & (experiments==expt)
        #if need to sample without replacement
        #sample_n = append(sample_n, sample(which(strat), as.integer(0.8*sum(strat)))
        sample_n = append(sample_n, sample(which(strat), sum(strat), replace=TRUE))
      }}
    for (clus in levels(clusters)){
      strat = (clusters==clus)
      #if need to sample without replacement
      #sample_p = append(sample_p, sample(which(strat), as.integer(0.8*sum(strat)))
      sample_p = append(sample_p, sample(which(strat), sum(strat), replace=TRUE))
      }
      slingshot_results_n <- getSlingshotResults(cellData[sample_n, ], conditions[sample_n], clusters[sample_p], start.cluster, end.clusters, pseudo.cond, pseudo.data[sample_p, ])}

    sub.pseudotimes_n <- slingshot_results_n$pseudotimes
    sub.weights_n <- slingshot_results_n$weights
    return(list(sample=sample_n, pseudotimes=sub.pseudotimes_n, weights = sub.weights_n))
  }

  if (parallel){
    BPPARAM = BiocParallel::bpparam()
    if (!is.null(n_cores)){
      BPPARAM$workers=as.integer(n_cores)}
    sampling_results = BiocParallel::bplapply(seq_len(nSamples), getSampling, cellData=cellData, conditions=conditions, clusters=clusters, start.cluster=start.cluster, end.clusters=end.clusters, pseudo.cond=pseudo.cond, pseudo.data = pseudo.data, BPPARAM = BPPARAM)
  }
  else {sampling_results = lapply(seq_len(nSamples), getSampling, cellData=cellData, conditions=conditions, clusters=clusters, start.cluster=start.cluster, end.clusters=end.clusters, pseudo.cond=pseudo.cond, pseudo.data = pseudo.data)}


  for (n in 1:length(sampling_results)){
    cur_sample = sampling_results[[n]]
    samples[, n] = cur_sample$sample
    sub.pseudotimes[, , n] = cur_sample$pseudotimes
    sub.weights[, , n] = cur_sample$weights
  }

  return(list(ori.pseudotimes=ori.pseudotimes, ori.weights=ori.weights, samples=samples, sub.pseudotimes=sub.pseudotimes, sub.weights=sub.weights))
}
