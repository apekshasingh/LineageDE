#' chooseK
#'
#' Fits alternative NB-GAM model for a random sample (n=100 by default) of genes
#' over different values of k (number of knots) and returns average AIC.  Minimum
#' AIC suggests value of k for downstream analysis.
#'
#' @param counts matrix (size = cells by genes) of raw (integer) count values of gene expression
#' @param pseudotimes matrix (size = cells by lineages) of pseudotime values
#' @param conditions factor (length = cells) specifying condition to which each cell belongs
#' @param weights matrix (size = cells by lineages) of weights from pseudotime analysis if applicable (otherwise NULL)
#' @param experiments factor (length = cells) specifying replicate/batch to which each cell belongs if applicable (otherwise NULL)
#' @param linCompare integer specifying which lineage (index) to test (defaults to 1 - only one lineage provided in pseudotimes)
#' @param n integer number of genes to randomly sample (default = 100)
#' @param k_vals vector of integers specifying number of knots in NB-GAM fit to trial (default = c(1:9))
#' @param genes vector of gene names or indices to sample from if desired (defaults to NULL and samples from all genes)
#' @param parallel boolean specifying whether or not to run in parallel (requires BiocParallel)
#' @param n_cores integer specifying number of cores to use if running in parallel
#' @param save.name string specifying file name if image of AIC results to be saved (if set to NULL, no image saved). Defaults to "k_results.png".
#' @param fit_time time permitted (in seconds) for NB-GAM fits for single gene (defaults to 100)
#' @param seed set random seed
#'
#' @returns a list containing AIC, R2 (r squared), and deviance explained values (matrix size = n by length(k_vals))
#'
#' @export
chooseK<- function(counts, pseudotimes, conditions, weights=NULL, experiments=NULL, linCompare=1, n=100, k_vals = 1:9, genes=NULL, parallel=FALSE, n_cores=NULL, save.name="k_results.png", fit_time=100, seed=123){

  ## DATA INPUTS
  # counts: matrix (cell by genes - raw values)

  # pseudotime results (no sampling input required here, just original):
  # pseudotimes = cells x lineages
  # weights = cells x lineages
  # if no weights provided, assume all cells have equal weight
  if (is.null(weights)){weights = array(1, dim(pseudotimes))}
  # all pseudotime samples should be normalized! (max 1) done in code below if needed

  #conditions: factor specifying the condition to which each cell belongs
  conditions = droplevels(conditions)

  #experiments: factor specifying the experiment/batch to which each cell belongs
  if (!is.null(experiments)){experiments = droplevels(experiments)}

  #max-normalize pseudotimes
  if (max(pseudotimes, na.rm=TRUE)>1){pseudotimes = pseudotimes/max(pseudotimes, na.rm=TRUE)
  warning("pseudotimes were not normalized, dividing by maximum value")}

  # linCompare specifies the index of the lineage being tested
  pseudotimes = pseudotimes[, linCompare]
  weights = weights[, linCompare]

  ## Specifications for chooseK
  # n: number of genes to sample to choose K
  # k_vals: range of k values (# knots in nbgam) to explore
  # genes: names/indices of genes to choose from for testing k val
  # parallel: TRUE or FALSE whether to fit genes in parallel (utilzing BiocParallel)
  # n_cores: Number of cores to utilize if running in parallel
  if (parallel){if(!is.numeric(n_cores)){stop("if setting parallel to true need to specify n_cores")}
    else if(n_cores%%1>0){stop("if setting parallel to true need to specify integer value for n_cores")}}

  # save.name: filepath specifying name of output (if NULL no image generated)
  # fit_time: time permitted (sec) for each gene model fitting (can increase if many errors)
  # seed: set for random sampling so genes selected is reproducible
  set.seed(seed)

  ## setting up additional input for fitting nbgam models

  ## Offset
  # will need to determine offset for each cell = log(Ni), where Ni is the libsize
  libSize = Matrix::rowSums(counts)

  select = libSize>0 #should have all cells with at least one count
  libSize = libSize[select]
  counts = counts[select, ]
  pseudotimes = pseudotimes[select]
  weights = weights[select]
  conditions = conditions[select]
  experiments = experiments[select]

  offsets = log(libSize)

  fit_nbgam<- function(gene_counts, pseudotimes, weights, conditions, experiments, offsets, k){

    min.pseudotimes = sapply(levels(conditions), function(condition){min(pseudotimes[conditions==condition & weights > 0], na.rm=TRUE)})
    commonStart = max(min.pseudotimes) #gives first knot location, ensures reasonable comparison between starts of lineages of differing lengths
    max.pseudotimes = sapply(levels(conditions), function(condition){max(pseudotimes[conditions==condition & weights > 0], na.rm=TRUE)})
    commonEnd = min(max.pseudotimes) #gives last knot location, ensures reasonable comparison between ends of lineages of differing lengths;
    #if final knot were after end of shorter lineage--would not be determined from data
    knot_locs = stats::quantile(pseudotimes[weights > 0 & pseudotimes<=commonEnd & pseudotimes>=commonStart], probs = (0:(k-1))/(k-1))
    if (length(unique(knot_locs))!=k){knot_locs = seq(commonStart, commonEnd, length=k)}
    knotList <- list(pseudotimes=knot_locs)

    ## FIT NBGAM w/ mgcv

    #fit alternative model (lineage specific coefficients)
    if(is.null(experiments)){if (k>2){smoothForm <- stats::as.formula("gene_counts ~ conditions + s(pseudotimes, by=conditions, bs='cr', k=k, fx=TRUE) + offset(offsets)")}
      else if(k==2){smoothForm <- stats::as.formula("gene_counts ~ conditions + conditions:pseudotimes + offset(offsets)")}
      else{smoothForm <- stats::as.formula("gene_counts ~ conditions + offset(offsets)")}}
    else{if (k>2){smoothForm <- stats::as.formula("gene_counts ~ conditions + experiments + s(pseudotimes, by=conditions, bs='cr', k=k, fx=TRUE) + offset(offsets)")}
      else if(k==2){smoothForm <- stats::as.formula("gene_counts ~ conditions + experiments + conditions:pseudotimes+ offset(offsets)")}
      else{smoothForm <- stats::as.formula("gene_counts ~ conditions + experiments + offset(offsets)")}}
    m <- try(mgcv::gam(smoothForm, family = "nb", knots = knotList, weights = weights/mean(weights[weights > 0]), control = mgcv::gam.control()))
    if (class(m)[1]=="try-error"){
      print("error fitting a model")
      m=NULL}
    return(m)}

  fit_nbgam_timed<-function(gene_counts, pseudotimes, weights, conditions, experiments, offsets, k, time){
    m<-try(R.utils::withTimeout(fit_nbgam(gene_counts, pseudotimes, weights, conditions, experiments, offsets, k), timeout=time, cpu=Inf, elapsed=time))
    if (class(m)[1]=="try-error"){
      print("error fitting a model - timeout")
      m=NULL}
    return(m)}

  # sample genes for testing k values
  if (is.null(genes)){genes = colnames(counts)}
  genes = sample(genes, n, replace=FALSE)

  evalK<-function(g, counts, pseudotimes, weights, conditions, experiments, offsets){
    # fit nbgam for range of k_vals; return fit model deviance, aic
    print(g)
    #will only keep cells that have some assignment to the lineage being compared (weight >  0)
    keep_cell = (!is.na(pseudotimes) & !is.na(weights)) & (weights > 0)
    gene_counts = counts[keep_cell, g]
    r2 = c()
    dev_exp = c()
    deviances = c()
    aic = c()

    for (k in k_vals){
      model = fit_nbgam_timed(gene_counts, pseudotimes[keep_cell], weights[keep_cell], conditions[keep_cell], experiments[keep_cell], offsets[keep_cell], k, fit_time)
      if(is.null(model)){r2 = c(r2, rep(NA, length(k_vals)-which(k==k_vals)+1))
        dev_exp = c(dev_exp, rep(NA, length(k_vals)-which(k==k_vals)+1))
        aic = c(aic, rep(NA, length(k_vals)-which(k==k_vals)+1))
        deviances = c(deviances, rep(NA, length(k_vals)-which(k==k_vals)+1))
        break}
      smry = summary(model)
      r2[length(r2)+1] = smry$r.sq
      dev_exp[length(dev_exp)+1] = smry$dev.expl

      deviances[length(deviances)+1] = model$deviance
      aic[length(aic)+1] = stats::AIC(model)}

    return(c(r2, dev_exp, deviances, aic))}

  if (parallel){
    BPPARAM = BiocParallel::bpparam()
    if (!is.null(n_cores)){
      BPPARAM$workers=as.integer(n_cores)}
    results = BiocParallel::bplapply(genes, evalK, counts=counts, pseudotimes=pseudotimes, weights=weights, conditions=conditions, experiments=experiments, offsets=offsets, BPPARAM = BPPARAM)}
  else{results = lapply(genes, evalK, counts=counts, pseudotimes=pseudotimes, weights=weights, conditions=conditions, experiments=experiments, offsets=offsets)}

  results = matrix(unlist(results), nrow=length(genes), ncol=4*length(k_vals), byrow=TRUE)
  rownames(results) = genes
  colnames(results) = rep(k_vals, 4)

  R2 = results[, 1:length(k_vals)]
  DevExpl = results[, (1+length(k_vals)):(2*length(k_vals))]
  Deviance = results[, (1+(2*length(k_vals))):(3*length(k_vals))]
  AIC = results[, (1+(3*length(k_vals))):(4*length(k_vals))]
  results = list(AIC=AIC, R2=R2, DevExpl=DevExpl, Deviance=Deviance)
  if (!is.null(save.name)){
    grDevices::png(file=save.name)
    plot_AIC = AIC/AIC[, 1]
    plot(k_vals, apply(plot_AIC, 2, mean), pch=19, ylab="AIC", xlab="k")
    graphics::lines(k_vals, apply(plot_AIC, 2, mean))
    graphics::arrows(k_vals, apply(plot_AIC, 2, mean) - apply(AIC, 2, stats::sd), apply(AIC, 2, mean) + apply(AIC, 2, stats::sd), length=0.05, angle=90, code=3)
    grDevices::dev.off()}

  print("choosing K done!")
  return(results)
}
