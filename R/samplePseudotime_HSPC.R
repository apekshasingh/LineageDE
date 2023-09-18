#' samplePseudotime_HSPC
#'
#' modified version of samplePseudotime to break derived pseudotime trajectories for HSPCs
#' into a HSC and primed progenitor lineages.  Furthermore adjust weights to account for differences
#' in projection distances to derived trajectories.
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
samplePseudotime_HSPC<- function(cellData, conditions, clusters, start.cluster, end.clusters, experiments, nSamples=100, pseudo.cond=NULL, pseudo.data = NULL, parallel=FALSE, n_cores=NULL, seed=123){

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
      br_pt<-get_br_pt(slingshot_results, cellData, nLin)
    }
    else if(!is.null(pseudo.cond)){
      ref_cells = conditions==pseudo.cond
      #obtain pseudotime reference
      slingshot_results<-slingshot::slingshot(cellData[ref_cells, ], clusters[ref_cells], start.clus = start.cluster, end.clus=end.clusters)
      max_pt = max(slingshot::slingPseudotime(slingshot_results), na.rm=TRUE)
      br_pt<-get_br_pt(slingshot_results, cellData[ref_cells, ], nLin)
      #project all cells onto reference
      slingshot_results<-slingshot::predict(slingshot_results, cellData)
    }
    else{
      #obtain pseudotime reference (alternative dataset)
      slingshot_results<-slingshot::slingshot(pseudo.data, clusters, start.clus = start.cluster, end.clus=end.clusters)
      max_pt = max(slingshot::slingPseudotime(slingshot_results), na.rm=TRUE)
      br_pt<-get_br_pt(slingshot_results, cellData[ref_cells, ], nLin)
      #project cells onto reference
      slingshot_results<-slingshot::predict(slingshot_results, rbind(cellData, pseudo.data))
    }

    #lineages = slingshot_results@lineages
    #curves = slingshot_results@curves}

    lineages = slingshot_results@metadata$lineages
    curves = slingshot_results@metadata$curves

    for (l in 1:nLin){
      lin_id = utils::tail(lineages[[l]], 1)
      lin_id = which(lin_id == end.clusters)

      pseudotimes[, lin_id] = curves[[l]]$lambda[1:nCell]
      weights[, lin_id] = curves[[l]]$w[1:nCell]
      distances[, lin_id] = curves[[l]]$dist_ind[1:nCell]}

    #max-normalize pseudotimes
    pseudotimes = pseudotimes/max_pt

    weights_results<-getNormWeights(pseudotimes, weights, distances, br_pt)
    pseudotimes <- weights_results$pseudotimes
    weights <- weights_results$weights

    return(list(pseudotimes=pseudotimes, weights=weights))
  }

  interp_curve <- function(time, curve){
    result = matrix(nrow=dim(curve)[2], ncol=1001)
    for (i in 1:dim(curve)[2]){
      sp = stats::smooth.spline(time, curve[, i])
      result[i, ] = stats::predict(sp, seq(from=0, to=1, by=0.001))$y}
    return(result)}

  max_traj_dist<-function(times, interp_curves){
    inter_traj_dist = c()
    for (i in 1:dim(interp_curves)[3]){
      inter_traj_dist = append(inter_traj_dist, max(stats::dist(interp_curves[, , i], method="euclidean")))}
    return(inter_traj_dist)}

  deriv_traj_dist<-function(times, inter_traj_dist){
    delta = times[2]-times[1]
    npts = length(inter_traj_dist)
    cfd = (inter_traj_dist[-c(1, 2)]-inter_traj_dist[-c(npts, npts-1)])/(2 * delta)
    return(cfd)}

  get_br_pt<-function(slingshot_results, cellData, nLin){
      lineages = slingshot_results@metadata$lineages
      curves = slingshot_results@metadata$curves
      max_pt = max(slingshot::slingPseudotime(slingshot_results), na.rm=TRUE)

      interp_curves = array(dim=c(nLin, ncol(cellData), 1001))
      for (l in 1:nLin){
        lin_id = utils::tail(lineages[[l]], 1)
        lin_id = which(lin_id == end.clusters)

        #obtain time associated with each curve
        curve_times = slingshot::predict(slingshot_results, rbind(cellData, curves[[l]]$s))
        curve_times = curve_times@metadata$curves[[l]]$lambda[-seq_len(nrow(cellData))]
        interp_curves[lin_id, , ] = interp_curve(curve_times/max_pt, curves[[l]]$s)}

      times = seq(from=0, to=1, by=0.001)
      inter_traj_dist = max_traj_dist(times, interp_curves)
      cfd = deriv_traj_dist(times, inter_traj_dist)

      times = times[-c(1, length(times))]
      br_pt = times[which.max(cfd)]
      return(br_pt)}

  # Scale Weights by Distances & Generate HSC vs Primed Lineages
  getNormWeights<- function(pseudotimes, weights, distances, br_pt, seed=123){

    set.seed(seed)

    dist_norm_wts = weights/sqrt(distances)
    wt = data.frame(dist_norm_wts)
    pt = data.frame(pseudotimes)
    colnames(wt)<-c("ErP", "MkP", "MyP", "LyP")
    colnames(pt)<-c("ErP", "MkP", "MyP", "LyP")

    #collapse the ErP and MkP lineages by taking higher wt
    EMP_pt = cbind(pt$ErP, pt$MkP)
    EMP_wt = cbind(wt$ErP, wt$MkP)
    EMP_select = max.col(EMP_wt)
    EMP_wt = EMP_wt[cbind(seq_len(dim(EMP_wt)[1]), EMP_select)]
    EMP_pt = EMP_pt[cbind(seq_len(dim(EMP_pt)[1]), EMP_select)]
    EMP_weights = EMP_wt
    EMP_weights[EMP_pt<=br_pt] = 0

    MyP_pt = pt$MyP
    MyP_wt = wt$MyP
    MyP_weights = MyP_wt
    MyP_weights[MyP_pt<=br_pt] = 0

    LyP_pt = pt$LyP
    LyP_wt = wt$LyP
    LyP_weights = LyP_wt
    LyP_weights[LyP_pt<=br_pt] = 0

    HSC_pt = cbind(EMP_pt, MyP_pt, LyP_pt)
    HSC_wt = cbind(EMP_wt, MyP_wt, LyP_wt)
    EMP_pt[EMP_pt<=br_pt] = NA
    MyP_pt[MyP_pt<=br_pt] = NA
    LyP_pt[LyP_pt<=br_pt] = NA

    HSC_select = HSC_pt<=br_pt
    num_HSC_select = rowSums(HSC_select, na.rm=TRUE)

    HSC_pt = rowSums(HSC_select*HSC_pt, na.rm=TRUE)/num_HSC_select
    HSC_wt = rowSums(HSC_select*HSC_wt, na.rm=TRUE)/num_HSC_select

    HSC_wt[is.na(HSC_wt)] = 0
    HSC_pt[is.na(HSC_pt)] = NA

    pseudotimes = cbind(HSC_pt, EMP_pt, MyP_pt, LyP_pt)
    weights = cbind(HSC_wt, EMP_weights, MyP_weights, LyP_weights)

    return(list(pseudotimes = pseudotimes, weights = weights))}

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
    else{
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
    return(list(sample=sample_n, pseudotimes=sub.pseudotimes_n, weights = sub.weights_n))}

  if (parallel){
    BPPARAM = BiocParallel::bpparam()
    if (!is.null(n_cores)){
      BPPARAM$workers=as.integer(n_cores)}
    sampling_results = BiocParallel::bplapply(seq_len(nSamples), getSampling, cellData=cellData, conditions=conditions, clusters=clusters, start.cluster=start.cluster, end.clusters=end.clusters, pseudo.cond=pseudo.cond, pseudo.data=pseudo.data, BPPARAM = BPPARAM)
  }
  else {sampling_results = lapply(seq_len(nSamples), getSampling, cellData=cellData, conditions=conditions, clusters=clusters, start.cluster=start.cluster, end.clusters=end.clusters, pseudo.cond=pseudo.cond, pseudo.data=pseudo.data)}

  for (n in 1:length(sampling_results)){
    cur_sample = sampling_results[[n]]
    samples[, n] = cur_sample$sample
    sub.pseudotimes[, , n] = cur_sample$pseudotimes
    sub.weights[, , n] = cur_sample$weights}

  return(list(ori.pseudotimes=ori.pseudotimes, ori.weights=ori.weights, samples=samples, sub.pseudotimes=sub.pseudotimes, sub.weights=sub.weights))}
