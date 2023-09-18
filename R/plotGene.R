#' plotGene
#'
#' Plots expression for gene of interest over pseudotime for a single lineage.
#' Plots both CPM and NB-GAM fits separated by specified conditions.
#'
#' @param gene name or index of gene to plot
#' @param counts matrix (size = cells by genes) of raw (integer) count values of gene expression
#' @param pseudotimes matrix (size = cells by lineages) of pseudotime values
#' @param conditions factor (length = cells) specifying condition to which each cell belongs
#' @param weights matrix (size = cells by lineages) of weights from pseudotime analysis if applicable (otherwise NULL)
#' @param experiments factor (length = cells) specifying replicate/batch to which each cell belongs if applicable (otherwise NULL)
#' @param linCompare integer specifying which lineage (index) to plot (defaults to 1 - only one lineage provided in pseudotimes)
#' @param k integer specifying number of knots in NB-GAM fit (default = 3)
#' @param seed set random seed
#'
#' @returns ggplot object
#'
#' @export
plotGene<-function(gene, counts, pseudotimes, conditions, weights=NULL, experiments=NULL, linCompare=1, k=3, seed=123){

  ## DATA INPUTS:

  # gene: names/index of gene to plot
  # counts: matrix (cell by genes - raw values)

  # pseudotime results:
  # pseudotimes = cells x lineages
  # weights = cells x lineages
  # if no weights provided, assume all cells have equal weight
  if (is.null(weights)){weights = array(1, dim(pseudotimes))}
  # all pseudotime samples should be normalized! (max 1)

  # conditions: factor specifying the condition of each cell
  conditions = droplevels(conditions)

  #experiments: factor specifying the experiment/batch to which each cell belongs
  if (!is.null(experiments)){experiments = droplevels(experiments)}

  # linCompare provides the index of the lineage to be plotted
  pseudotimes = pseudotimes[, linCompare]
  weights = weights[, linCompare]

  # k is the number of knots used in the nbgam
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
  min.pseudotimes = sapply(levels(conditions), function(condition){min(pseudotimes[conditions==condition & weights > 0], na.rm=TRUE)})
  commonStart = max(min.pseudotimes) #gives first knot location, ensures reasonable comparison between starts of lineages of differing lengths
  max.pseudotimes = sapply(levels(conditions), function(condition){max(pseudotimes[conditions==condition & weights > 0], na.rm=TRUE)})
  commonEnd = min(max.pseudotimes) #gives last knot location, ensures reasonable comparison between ends of lineages of differing lengths;
  # if final knot were after end of shorter lineage--would not be determined from data
  knot.locs = stats::quantile(pseudotimes[weights > 0 & pseudotimes<=commonEnd & pseudotimes>=commonStart], probs = (0:(k-1))/(k-1))
  if (length(unique(knot.locs))!=k){knot_locs = seq(commonStart, commonEnd, length=k)}
  knotList <- list(pseudotimes=knot.locs)

  getGeneFit<- function(gene_counts, pseudotimes, weights, conditions, experiments, offsets, knotList, k){
    ## fitting nbgam with mgcv, independent smooth per condition, penalized cubic regression
    if (is.null(experiments)){if(k>2){smoothForm <- stats::as.formula("gene_counts ~ conditions + s(pseudotimes, by=conditions, bs='cr', k=k, fx=TRUE) + offset(offsets)")}
      else if(k==2){smoothForm <- stats::as.formula("gene_counts ~ conditions + conditions:pseudotimes + offset(offsets)")}
      else{smoothForm <- stats::as.formula("gene_counts ~ conditions + offset(offsets)")}}
    else{if(k>2){smoothForm <- stats::as.formula("gene_counts ~ conditions + experiments + s(pseudotimes, by=conditions, bs='cr', k=k, fx=TRUE) + offset(offsets)")}
      else if(k==2){smoothForm <- stats::as.formula("gene_counts ~ conditions + experiments + conditions:pseudotimes + offset(offsets)")}
      else{smoothForm <- stats::as.formula("gene_counts ~ conditions + experiments + offset(offsets)")}}
    s <- mgcv::s
    m <- mgcv::gam(smoothForm, family = "nb", knots = knotList, weights = weights/mean(weights[weights >0]), control = mgcv::gam.control())
    return(m)}

  modelFit = getGeneFit(counts[,gene], pseudotimes, weights, conditions, experiments, offsets, knotList, k)

  CPM = (counts[,gene]/libSize)*1e6
  groups = factor(paste(experiments, conditions))
  data = data.frame(pseudotimes, weights, groups, CPM)
  colors = RColorBrewer::brewer.pal(length(levels(groups)), name="Set1")
  suppressWarnings(p <-ggplot2::ggplot(data=data, mapping=ggplot2::aes(x=pseudotimes, y=CPM, colour=groups))+ggplot2::geom_point(mapping=ggplot2::aes(size=weights, alpha=weights))+ggplot2::scale_size(range=c(0, 2))+ggplot2::scale_alpha(range=c(0, 0.1))+ggplot2::scale_color_manual(values =colors))

  if(is.null(experiments)){
  for (cond in levels(conditions)){
    select = conditions==cond
    model_CPM = (exp(stats::predict(modelFit, data.frame(pseudotimes=pseudotimes[select], conditions = rep(as.character(cond), sum(select)), offsets=rep(0, sum(select))))))*1e6
    c = which(cond == levels(groups))
    data = data.frame(pseudotimes = pseudotimes[select], model_CPM = model_CPM)
    suppressWarnings(p <- p +ggplot2::geom_line(mapping=ggplot2::aes(x=pseudotimes, y=model_CPM), colour=colors[c], linewidth = 1, data=data))}}
  else{
  for (expt in levels(experiments)){
    for (cond in levels(conditions)){
      select = conditions==cond & experiments==expt
      model_CPM = (exp(stats::predict(modelFit, data.frame(pseudotimes=pseudotimes[select], conditions = rep(as.character(cond), sum(select)), experiments = rep(as.character(expt), sum(select)), offsets=rep(0, sum(select))))))*1e6
      c = which(paste(expt, cond) == levels(groups))
      data = data.frame(pseudotimes = pseudotimes[select], model_CPM = model_CPM)
      suppressWarnings(p <- p +ggplot2::geom_line(mapping=ggplot2::aes(x=pseudotimes, y=model_CPM), colour=colors[c], linewidth = 1, data=data))}}}

  return(suppressWarnings(print(p)))}
