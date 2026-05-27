#' The simPICcount class
#' 
#' S4 class that holds parameters for simPIC simulation.
#' 
#' @section Parameters:
#' simPIC simulation parameters:
#' \describe{
#'     \item{\code{nPeaks}}{The number of peaks to simulate.}
#'     \item{\code{nCells}}{The number of cells to simulate.}
#'     \item{\code{[seed]}}{Seed to use for generating random numbers.}
#'     \item{\code{[default]}}{Logical value indicating whether to use default
#'     parameters (TRUE) or learn parameters from data (FALSE).}
#'     \item{\emph{Library size parameters}}{
#'         \describe{
#'             \item{\code{lib.size.meanlog}}{meanlog (location) parameter
#'             for the library size log-normal distribution.}
#'             \item{\code{lib.size.sdlog}}{sdlog (scale) parameter for the
#'             library size log-normal distribution.}
#'  }
#' }
#'      \item{\emph{Peak mean parameters}}{
#'         \describe{
#'             \item{\code{mean.scale}}{scale parameter for the mean
#'              weibull distribution.}
#'             \item{\code{mean.shape}}{shape parameter for the mean
#'             weibull distribution.}
#'  }
#' }
#'      \item{\emph{Cell sparsity parameters}}{
#'         \describe{
#'             \item{\code{sparsity}}{Probability that contributes to the
#'             sparsity of the final simulated matrix.}
#'  }
#' }
#' }
#'
#' @return a simPIC class object.
#' The parameters not shown in brackets can be estimated from real data
#' using \code{\link{simPICestimate}}. For details of the simPIC simulation
#' see \code{\link{simPICsimulate}}. The default parameters are based on the
#' PBMC10k dataset and can be reproduced using the test data and script
#' provided in \code{inst/scripts}.
#' @name simPICcount
#' @rdname simPICcount
#' @exportClass simPICcount
#' @aliases simPICcount-class
setClass("simPICcount",
    slots = c(
        nPeaks = "numeric",
        nCells = "numeric",
        seed = "numeric",
        default = "logical",
        pm.distr = "character",
        lib.size.meanlog = "numeric",
        lib.size.sdlog = "numeric",
        peak.mean.shape = "numeric",
        peak.mean.rate = "numeric",
        peak.mean.scale = "numeric",
        peak.mean.pi = "numeric",
        peak.mean.meanlog = "numeric",
        peak.mean.sdlog = "numeric",
        sparsity = "numeric",
        nGroups = "numeric",
        group.prob = "numeric",
        da.prob = "numeric",
        da.downProb = "numeric",
        da.facLoc = "numeric",
        da.facScale = "numeric",
        bcv.common = "numeric",
        bcv.df = "numeric",
        nBatches = "numeric",
        batchCells = "numeric",
        batch.facLoc = "numeric",
        batch.facScale = "numeric",
        batch.rmEffect = "logical"
    ),
    prototype = prototype(
        nPeaks = 5000,
        nCells = 700,
        seed = sample(seq_len(1e5), 1),
        default = TRUE,
        pm.distr = "lngamma",
        lib.size.meanlog = 6.687082,
        lib.size.sdlog = 0.344361,
        peak.mean.shape = 0.7909301,
        peak.mean.rate = 7.100648,
        peak.mean.scale = 0.09522228,
        peak.mean.pi = -17.17441,
        peak.mean.meanlog = -2.825233,
        peak.mean.sdlog = -1.366378,
        sparsity = replicate(200000, runif(1, min = 0, max = 0.99)),
        nGroups = 1,
        group.prob = 1,
        da.prob = 0.13,
        da.downProb = 0.5,
        da.facLoc = 0.1,
        da.facScale = 0.4,
        bcv.common = 0.1,
        bcv.df = 60,
        nBatches = 1,
        batchCells = 700,
        batch.facLoc = 0.1,
        batch.facScale = 0.1,
        batch.rmEffect = FALSE
    )
)
