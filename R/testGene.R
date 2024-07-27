#' testGene
#'
#' Fits NB-GAM to test alternative vs null model for a specific gene.  Helper function for LineageDE.
#' Performs likelihood ratio test and returns p-value based off of chi-squared distribution (ori.pVal)
#' If samples of pseudotime are provided, will also perform permutation test and give an empirical p-value (emp.pVal)
#' and a parametric p-value (par.pVal) from fitting null distribution with a gamma distribution
#'
#' @param gene names or index specifying gene to be tested
#' @param counts matrix (size = cells by genes) of raw (integer) count values of gene expression
#' @param ori.pseudotimes matrix (size = cells by lineages) of pseudotime values
#' @param conditions factor (length = cells) specifying condition to which each cell belongs
#' @param ori.weights matrix (size = cells by lineages) of weights from pseudotime analysis if applicable (otherwise NULL)
#' @param experiments factor (length = cells) specifying replicate/batch to which each cell belongs if applicable (otherwise NULL)
#' @param offsets offset used in fitting NB-GAM based off of library size for each cell
#' @param knotList list specifying knot locations used in NB-GAM fit
#' @param k integer specifying number of knots in NB-GAM fit
#' @param samples matrix (size = cells by num_samples) specifying indices cells utilized in each sample of pseudotime if applicable (otherwise NULL)
#' @param sub.pseudotimes matrix (size = cells by lineages by num_samples) samples of pseudotime values if applicable (otherwise NULL)
#' @param sub.weights matrix (size = cells by lineages by num_samples) samples of pseudotime weights if applicable (otherwise NULL)
#' @param fit_time time permitted (in seconds) for NB-GAM fit for single gene
#' @param seed set random seed
#'
#' @returns LRT, ori.pVal, (emp.pVal), (par.pVal) for gene tested
#'
#' @export
testGene<- function(gene, counts, ori.pseudotimes, ori.weights, conditions, experiments, offsets, knotList, k, samples, sub.pseudotimes, sub.weights, fit_time, seed){

  set.seed(seed)

  ## NBGAM MODEL FITTING
  fit_nbgam<- function(gene_counts, pseudotimes, weights, conditions, experiments, offsets, knotList, k, fit_time){

    ## fitting nbgam with mgcv

    #fit null model (common coefficients across lineages)
    if (is.null(experiments)){if(k>2){smoothForm <- stats::as.formula("gene_counts ~ s(pseudotimes, bs='cr', k=k, fx=TRUE) + offset(offsets)")}
      else if(k==2){smoothForm <- stats::as.formula("gene_counts ~ pseudotimes + offset(offsets)")}
      else{smoothForm <- stats::as.formula("gene_counts ~ 1 + offset(offsets)")}}
    else{if(k>2){smoothForm <- stats::as.formula("gene_counts ~ experiments + s(pseudotimes, bs='cr', k=k, fx=TRUE) + offset(offsets)")}
      else if(k==2){smoothForm <- stats::as.formula("gene_counts ~ experiments + pseudotimes + offset(offsets)")}
      else{smoothForm <- stats::as.formula("gene_counts ~ experiments + offset(offsets)")}}
    s <- mgcv::s
    m0 <- try(R.utils::withTimeout(mgcv::gam(smoothForm, family = "nb", knots = knotList, weights = weights/mean(weights[weights>0]), control = mgcv::gam.control()),
                                   timeout=fit_time, cpu=Inf, elapsed=fit_time))
    if (class(m0)[1]=="try-error"){
      print("error fitting a null model")
      m0=NULL}

    #fit alternative model (lineage specific coefficients)
    if (is.null(experiments)){if(k>2){smoothForm <- stats::as.formula("gene_counts ~ conditions + s(pseudotimes, by=conditions, bs='cr', k=k, fx=TRUE) + offset(offsets)")}
      else if(k==2){smoothForm <- stats::as.formula("gene_counts ~ conditions + conditions:pseudotimes + offset(offsets)")}
      else{smoothForm <- stats::as.formula("gene_counts ~ conditions + offset(offsets)")}}
    else{if(k>2){smoothForm <- stats::as.formula("gene_counts ~ conditions + experiments + s(pseudotimes, by=conditions, bs='cr', k=k, fx=TRUE) + offset(offsets)")}
      else if(k==2){smoothForm <- stats::as.formula("gene_counts ~ conditions + experiments + conditions:pseudotimes + offset(offsets)")}
      else{smoothForm <- stats::as.formula("gene_counts ~ conditions + experiments + offset(offsets)")}}
    s <- mgcv::s
    m <- try(R.utils::withTimeout(mgcv::gam(smoothForm, family = "nb", knots = knotList, weights = weights/mean(weights[weights>0]), control = mgcv::gam.control()),
                                  timeout=fit_time, cpu=Inf, elapsed=fit_time))
    if (class(m)[1]=="try-error"){
      print("error fitting a alt model")
      m=NULL}
    model_fit<-list(null_model=m0, alt_model=m)
    return(model_fit)}

  ## FIT ORIGINAL DATA
  print(gene)
  gene_counts = counts[, gene]

  #will only keep cells that have some assignment to the lineage being compared (weight >= 0)
  keep_cell = (!is.na(ori.pseudotimes) & !is.na(ori.weights)) & (ori.weights>0) & (offsets>-Inf)

  ori.model_fit = fit_nbgam(gene_counts[keep_cell], ori.pseudotimes[keep_cell], ori.weights[keep_cell], conditions[keep_cell], experiments[keep_cell], offsets[keep_cell], knotList, k, fit_time)
  ori.LRT = try(-2*(stats::logLik(ori.model_fit$null_model)-stats::logLik(ori.model_fit$alt_model)), silent = TRUE)
  if (class(ori.LRT)=="try-error") {
    print("original model fitting failed or execeeded permitted time")
    ori.LRT=NA}
  ori.pVal = 1-stats::pchisq(ori.LRT, df = k*(length(levels(conditions))-1))

  # PERMUTE AND FIT SUBSAMPLED DATA
  getSampleLRT<-function(s, gene_counts, sub.pseudotimes, sub.weights, conditions, experiments, offsets, knotList, k){
    gene_counts.sample = gene_counts[samples[, s]]
    offsets.sample = offsets[samples[,s]]
    pseudotimes.sample = sub.pseudotimes[, s]
    weights.sample = sub.weights[, s]
    experiments.sample = experiments[samples[, s]]
    # for the permutation test we will permute the condition labels
    conditions.sample = sample(conditions[samples[, s]], length(conditions[samples[, s]]), replace=FALSE)
    #again will only keep cells that have some assignment to the lineage being compared (weight > 0)
    keep_cell = (!is.na(pseudotimes.sample) & !is.na(weights.sample)) & (weights.sample>0) & (offsets.sample>-Inf)

    sub.model_fit = fit_nbgam(gene_counts.sample[keep_cell], pseudotimes.sample[keep_cell], weights.sample[keep_cell], conditions.sample[keep_cell], experiments.sample[keep_cell], offsets.sample[keep_cell], knotList, k, fit_time)
    null.LRT = try(-2*(stats::logLik(sub.model_fit$null_model)-stats::logLik(sub.model_fit$alt_model)), silent = TRUE)
    if (class(null.LRT)=="try-error") {
      print("sample model fitting failed or execeeded permitted time")
      null.LRT=NA}
    return(null.LRT)}

  if (!is.null(samples) & !is.na(ori.LRT)){

    null.LRT = c()
    for (s in 1:ncol(samples)){
      null.LRT.sample = getSampleLRT(s, gene_counts, sub.pseudotimes, sub.weights, conditions, experiments, offsets, knotList, k)
      if (is.na(null.LRT.sample)){break}
      null.LRT = append(null.LRT, null.LRT.sample)
    }


    # calculate p-val from empirical null distribution
    emp.pVal = (sum((null.LRT>=ori.LRT), na.rm=TRUE)+1)/(sum(!is.na(null.LRT))+1)

    # calculate p-val from parametric estimate of empirical null distribution
    # fit gamma distribution using MLE and fitdistrplus package
    # dist_null.LRT <-try(fitdistrplus::fitdist(null.LRT[!is.na(null.LRT)], "gamma"), silent=TRUE)
    dist_null.LRT <-try(fitdistrplus::fitdist(null.LRT[(!is.na(null.LRT))&(null.LRT>0)], "gamma", lower=c(0, 0)), silent=TRUE)
    if (class(dist_null.LRT)=="try-error"){dist_null.LRT=NULL}
    shape = dist_null.LRT$estimate["shape"]
    rate = dist_null.LRT$estimate["rate"]

    par.pVal = try(1-stats::pgamma(ori.LRT, shape, rate), silent=TRUE)
    if (class(par.pVal)=="try-error"){par.pVal=NA}}

  else{emp.pVal = NA
  par.pVal = NA}


  gene_result = c(gene, ori.LRT, ori.pVal, emp.pVal, par.pVal)
  return(gene_result)}
