#' Example Input Data for LineageDE
#'
#' A list containing example input data to run LineageDE
#'
#' \describe{
#'   \item{counts}{matrix of gene expression counts (cells x genes)}
#'   \item{conditions}{factor specifying experimental conditions}
#'   \item{experiments}{factor specifying experimental replicates}
#'   \item{ori.pseudotimes}{matrix of pseudotime values (cells x lineages)}
#'   \item{ori.weights}{matrix of pseudotime weights (cells x lineages)}
#'   \item{samples}{matrix of cell indices sampled (cells x number_samples=100)}
#'   \item{sub.pseudotimes}{matrix of sampled pseudotime values (cells x lineages x 100)}
#'   \item{sub.weights}{matrix of sampled pseudotime weights (cells x lineages x 100)}
#'
#' }
#'
#' @docType data
#' @name ExampleInputData
#' @usage data(ExampleInputData)
NULL

#' Example Input for samplePsuedotime
#'
#' A list containing example input data to run samplePseudotime_HSPC
#'
#' \describe{
#'   \item{cellData}{matrix of PCA projection (cells x dimensions)}
#'   \item{conditions}{factor specifying experimental conditions}
#'   \item{experiments}{factor specifying experimental replicates}
#'   \item{clusters}{factor specifying input cluster for each cell (input for Slingshot)}
#' }
#'
#' @docType data
#' @keywords datasets
#' @name exPseudotimeInput
#' @usage data(exPseudotimeInput)
NULL

