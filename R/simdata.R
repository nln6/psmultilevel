#' Simulated multi-level data.
#'
#' This data set is simulated to demonstrate the use of the
#' package \code{psmultilevel}. The treatment is randomly assigned to 30\% of the
#' observations. The outcome is generated from a cluster-specific random effects
#' model. The true ATE/ ATT is fixed at -1.
#'
#' @docType data
#'
#' @usage data(simdata)
#'
#' @format \code{simdata} is a dataframe with 1000 observations (rows) and 5
#'   variables (columns) named \code{X1}, \code{X2}, \code{cluster}, \code{Z},
#'   and \code{Y}. The first two are independent observed variables.
#'   \code{cluster}  identifies the cluster to which the observation belongs.
#'   \code{Z} denotes the treatment assignment (1 = Treated, 0 = Untreated).
#'   \code{Y} is the outcome of interest.
#'
#' @keywords datasets
#'
#' @examples
#' data(simdata)
#' head(simdata)
#'
"simdata"
