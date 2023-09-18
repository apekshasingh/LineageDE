#' LineageDE
#'
#' Fits NB-GAM to test alternative vs null model for collection of genes.
#' Performs likelihood ratio test and returns p-value based off of chi-squared distribution (ori.pVal)
#' If samples of pseudotime are provided, will also perform permutation test and give an empirical p-value (emp.pVal)
#' and a parametric p-value (par.pVal) from fitting null distribution with a gamma distribution
#'
#' @param genes vector of names or indices specifying genes to be tested
#' @param counts matrix (size = cells by genes) of raw (integer) count values of gene expression
#' @param ori.pseudotimes matrix (size = cells by lineages) of pseudotime values
#' @param conditions factor (length = cells) specifying condition to which each cell belongs
#' @param ori.weights matrix (size = cells by lineages) of weights from pseudotime analysis if applicable (otherwise NULL)
#' @param experiments factor (length = cells) specifying replicate/batch to which each cell belongs if applicable (otherwise NULL)
#' @param linCompare integer specifying which lineage (index) to test (defaults to 1 - only one lineage provided in pseudotimes)
#' @param k integer specifying number of knots in NB-GAM fit (default = 3)
#' @param samples matrix (size = cells by num_samples) specifying indices cells utilized in each sample of pseudotime if applicable (otherwise NULL)
#' @param sub.pseudotimes matrix (size = cells by lineages by num_samples) samples of pseudotime values if applicable (otherwise NULL)
#' @param sub.weights matrix (size = cells by lineages by num_samples) samples of pseudotime weights if applicable (otherwise NULL)
#' @param parallel boolean specifying whether or not to run in parallel (requires BiocParallel)
#' @param n_cores integer specifying number of cores to use if running in parallel
#' @param output integer specifying how often to output results to csv file (number of genes).  Defaults to NULL (output after all genes completed).
#' @param filePath string specifying path of directory where files are to be saved (Defaults to current working directory)
#' @param save.name string specifying prefix name of saved csv files (Defaults to "DE_results"). Specify "None" if do not want to save output to csv.
#' @param fit_time time permitted (in seconds) for NB-GAM fits for single gene (defaults to 100)
#' @param seed set random seed
#'
#' @returns dataframe with LRT, ori.pVal, (emp.pVal), (par.pVal) for each gene tested
#'
#' @export
LineageDE<- function(genes, counts, ori.pseudotimes, conditions, ori.weights=NULL, experiments=NULL, linCompare=1, k=3, samples=NULL, sub.pseudotimes=NULL, sub.weights=NULL, parallel=FALSE, n_cores=NULL, output=NULL, filePath=NULL, save.name=NULL, fit_time=100, seed=123){

  ## DATA INPUTS:
  # assuming multiple subsamplings of pseudotime output provided as input
  # if empirical & parametric p-values are to be calculated

  # genes: names/indices of genes to be tested
  # counts: matrix (cell by genes - raw values)

  # pseudotime results:
  # ori.pseudotimes = cells x lineages
  # ori.weights = cells x lineages
  # samples = indices (80%-100% cells) x n_samples
  # sub.pseudotimes = cells (80%-100%) x lineages x n_samples
  # sub.weights = cells (80%-100%) x lineages x n_samples
  # if no weights provided, assume all cells have equal weight
  if (is.null(ori.weights)){ori.weights = array(1, dim(ori.pseudotimes))}
  if (is.null(sub.weights) & !is.null(samples)){sub.weights = array(1, dim(sub.pseudotimes))}
  # all pseudotime samples should be normalized! (max 1)

  # conditions: factor specifying the condition of each cell
  conditions = droplevels(conditions)

  #experiments: factor specifying the experiment/batch to which each cell belongs
  if (!is.null(experiments)){experiments = droplevels(experiments)}

  #max-normalize pseudotimes if not already performed
  if (max(ori.pseudotimes, na.rm=TRUE)>1){ori.pseudotimes = ori.pseudotimes/max(ori.pseudotimes, na.rm=TRUE)
  warning("ori.pseudotimes were not normalized, dividing by maximum value")}
  if (max(sub.pseudotimes, na.rm=TRUE)>1){sub.pseudotimes = sweep(sub.pseudotimes, 3, FUN="/", apply(sub.pseudotimes, 3, max, na.rm=TRUE))
  warning("sub.pseudotimes were not normalized, dividing by maximum value per sample")}

  # linCompare provides the index of the lineage to be tested
  ori.pseudotimes = ori.pseudotimes[, linCompare]
  ori.weights = ori.weights[, linCompare]
  if (!is.null(samples)){
  sub.pseudotimes = sub.pseudotimes[, linCompare, ]
  sub.weights = sub.weights[, linCompare, ]}

  # k is the number of knots used in the nbgam
  # filePath specifies where to save statistical testing results for each gene
  if (is.null(filePath)){filePath=getwd()}
  # save.name specifies name of csv file for results
  if (is.null(save.name)){save.name="DE_results"}
  # parallel: TRUE or FALSE whether to fit genes in parallel (utilzing BiocParallel)
  # n_cores: Number of cores to utilize if running in parallel
  if (parallel){if(!is.numeric(n_cores)){stop("if setting parallel to true need to specify n_cores")}
    else if(n_cores%%1>0){stop("if setting parallel to true need to specify integer value for n_cores")}}

  #output specifies how often to return results to csv (number of genes)
  if (is.null(output)){output=length(genes)}

  # fit_time: time permitted (sec) for each gene model fitting (can increase if many errors)
  # seed: set for any random sampling is reproducible
  set.seed(seed)

  ## setting up additional input for fitting nbgam models

  ## Offsets
  # will need to determine offset for each cell = log(Ni), where Ni is the libsize
  libSize = Matrix::rowSums(counts)
  offsets = log(libSize)
  gc()

  ## Determine knot locations
  # quantile based method
  min.pseudotimes = sapply(levels(conditions), function(condition){min(ori.pseudotimes[conditions==condition & ori.weights > 0], na.rm=TRUE)})
  commonStart = max(min.pseudotimes) #gives first knot location, ensures reasonable comparison between starts of lineages of differing lengths
  max.pseudotimes = sapply(levels(conditions), function(condition){max(ori.pseudotimes[conditions==condition & ori.weights > 0], na.rm=TRUE)})
  commonEnd = min(max.pseudotimes) #gives last knot location, ensures reasonable comparison between ends of lineages of differing lengths;
  # if final knot were after end of shorter lineage--would not be determined from data
  knot.locs = stats::quantile(ori.pseudotimes[ori.weights > 0 & ori.pseudotimes<=commonEnd & ori.pseudotimes>=commonStart], probs = (0:(k-1))/(k-1))
  if (length(unique(knot.locs))!=k){knot_locs = seq(commonStart, commonEnd, length=k)}
  knotList <- list(pseudotimes=knot.locs)

  geneSplits = split(genes, ceiling(seq_along(genes)/output))
  print("running LineageDE")

  results = data.frame(matrix(ncol=5, nrow=0))
  for (i in 1:length(geneSplits)){
    cur_genes = geneSplits[[i]]

    ## Parallelization Settings
    if (parallel){
      BPPARAM <- BiocParallel::bpparam()
      if (!is.null(n_cores)){
        BPPARAM$workers=as.integer(n_cores)}

      cur_results <- BiocParallel::bplapply(cur_genes, testGene,
                                            counts=counts,
                                            ori.pseudotimes=ori.pseudotimes,
                                            ori.weights=ori.weights,
                                            conditions=conditions,
                                            experiments=experiments,
                                            offsets=offsets,
                                            knotList=knotList,
                                            k=k,
                                            samples=samples,
                                            sub.pseudotimes=sub.pseudotimes,
                                            sub.weights=sub.weights,
                                            fit_time=fit_time,
                                            seed=seed,
                                            BPPARAM = BPPARAM)}

    else {cur_results <- lapply(cur_genes, testGene,
                                counts=counts,
                                ori.pseudotimes=ori.pseudotimes,
                                ori.weights=ori.weights,
                                conditions=conditions,
                                experiments=experiments,
                                offsets=offsets,
                                knotList=knotList,
                                k=k,
                                samples=samples,
                                sub.pseudotimes=sub.pseudotimes,
                                sub.weights=sub.weights,
                                fit_time=fit_time,
                                seed=seed)}

    cur_results = data.frame(matrix(unlist(cur_results), nrow=length(cur_genes), ncol=5, byrow=TRUE))
    colnames(cur_results) = c("gene", "LRT", "ori.pVal", "emp.pVal", "par.pVal")
    if (!(save.name=="None")){utils::write.csv(cur_results, file=paste0(filePath, "/", save.name, "_", i, ".csv"), row.names = FALSE)}
    results = rbind(results, cur_results)
    print(paste("finished gene split", i, "of", length(geneSplits)))}

  if (!is.null(samples)){results$pVal.adj = stats::p.adjust(results$par.pVal, method="BH")}
  else {results$pVal.adj = stats::p.adjust(results$ori.pVal, method="BH")}
  print("LineageDE done!")
  return(results)
}
