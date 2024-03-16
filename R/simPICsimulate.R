#' simPIC simulation
#'
#' Simulate peak by cell count matrix from a sparse single-cell ATAC-seq
#' peak by cell input using simPIC methods.
#'
#' @param object simPICcount object with simulation parameters.
#' See \code{\link{simPICcount}} for details.
#' @param pm.distr distribution parameter for peak means.
#' Available distributions: gamma, weibull, lngamma, pareto.
#' Default is weibull.
#' @param verbose logical variable. Prints the simulation progress if TRUE.
#' @param ... Any additional parameter settings to override what is provided
#'        in \code{simPICcount} object.
#'
#' @return SingleCellExperiment object containing the simulated counts.
#'
#' @details
#' simPIC provides the option to manually adjust each of the
#' \code{simPICcount} object parameters by calling
#' \code{\link{setsimPICparameters}}.
#'
#'  The simulation involves following steps:
#'  \enumerate{
#'  \item Set up simulation parameters
#'  \item Set up SingleCellExperiment object
#'  \item Simulate library sizes
#'  \item Simulate sparsity
#'  \item Simulate peak means
#'  \item Create final synthetic counts
#'  }
#'
#' The final output is a
#' \code{\link[SingleCellExperiment]{SingleCellExperiment}} object that
#' contains the simulated count matrix. The parameters are stored in the
#' \code{\link{colData}} (for cell specific information),
#' \code{\link{rowData}} (for peak specific information) or
#' \code{\link{assays}} (for peak by cell matrix) slots. This additional
#' information includes:
#' @examples
#' # default simulation
#' sim <- simPICsimulate(pm.distr = "weibull")
#'
#' @importFrom SummarizedExperiment rowData rowData<-
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom SummarizedExperiment assays assays<-
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @export
simPICsimulate <- function(object = newsimPICcount(),
                            verbose = TRUE, pm.distr = "weibull", ...) {
    checkmate::assertClass(object, "simPICcount")
    if (verbose) {
        message("simPIC is:")
        message("updating parameters...")
    }
    object <- setsimPICparameters(object, ...)
    validObject(object)

    seed <- simPICget(object, "seed")
    nCells <- simPICget(object, "nCells")
    nPeaks <- simPICget(object, "nPeaks")

    if (verbose) {
        message("setting up SingleCellExperiment object...")
    }
    cell.names <- paste0("Cell", seq_len(nCells))
    peak.names <- paste0("Peak", seq_len(nPeaks))
    cells <- data.frame(Cell = cell.names)
    rownames(cells) <- cell.names
    peaks <- data.frame(Peak = peak.names)
    rownames(peaks) <- peak.names

    sim <- SingleCellExperiment(
        rowData = peaks, colData = cells, metadata = list(Params = object)
    )
    if (verbose) {
        message("Simulating library size...")
        sim <- simPICsimulateLibSize(object, sim)
        gc()
        message("Simulating peak mean...")
        sim <- simPICsimulatePeakMean(object, sim, pm.distr)
        gc()
        message("Simulating true counts...")
    }
    sim <- simPICsimulateTrueCounts(object, sim)
    rownames(BiocGenerics::counts(sim)) <- peak.names
    colnames(BiocGenerics::counts(sim)) <- cell.names
    if (verbose) {
        message("Done!!")
    }
    return(sim)
}

#' Simulate simPIC library sizes
#'
#' Generate library sizes for cells in simPIC simulation based on the
#' estimated values of mus and sigmas.
#'
#' @param sim SingleCellExperiment object containing simulation parameters.
#' @param object simPICcount object with simulation parameters.
#' @param verbose logical. To print progress messages.
#'
#' @return SingleCellExperiment object with simulated library sizes.
#'
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom stats rlnorm
simPICsimulateLibSize <- function(object, sim, verbose) {
    nCells <- simPICget(object, "nCells")
    lib.size.meanlog <- simPICget(object, "lib.size.meanlog")
    lib.size.sdlog <- simPICget(object, "lib.size.sdlog")

    lib.size <- rlnorm(
        n = nCells,
        meanlog = lib.size.meanlog,
        sdlog = lib.size.sdlog
    )

    lib.size <- round(lib.size)
    colData(sim)$exp.libsize <- lib.size
    return(sim)
}

#' Simulate simPIC peak means.
#'
#' Generate peak means for cells in simPIC simulation based on the estimated
#' values of shape and rate parameters.
#'
#' @param sim SingleCellExperiment object containing simulation parameters.
#' @param object simPICcount object with simulation parameters.
#' @param pm.distr distribution parameter for peak means.
#' Available distributions: gamma, weibull, lngamma, pareto.
#' Default is weibull.
#' @param verbose logical. Whether to print progress messages.
#'
#' @return SingleCellExperiment object with simulated peak means.
#'
#' @importFrom SummarizedExperiment rowData rowData<-
#' @importFrom stats rgamma rweibull dgamma dlnorm pgamma plnorm runif
#' @importFrom Matrix sparseMatrix
#' @importFrom actuar rpareto
simPICsimulatePeakMean <- function(object, sim, pm.distr, verbose) {
    nPeaks <- simPICget(object, "nPeaks")
    message("using ", pm.distr, " distribution for simulating peak mean")

    switch(pm.distr,
        gamma = {
            peak.mean.rate <- simPICget(object, "peak.mean.rate")
            peak.mean.shape <- simPICget(object, "peak.mean.shape")
            peak.means <- rgamma(
                n = nPeaks, shape = peak.mean.shape,
                rate = peak.mean.rate
            )
        },
        weibull = {
            peak.mean.shape <- simPICget(object, "peak.mean.shape")
            peak.mean.scale <- simPICget(object, "peak.mean.scale")
            peak.means <- rweibull(
                n = nPeaks, shape = peak.mean.shape,
                scale = peak.mean.scale
            )
        },
        pareto = {
            peak.mean.shape <- simPICget(object, "peak.mean.shape")
            peak.mean.scale <- simPICget(object, "peak.mean.scale")
            peak.means <- rpareto(
                n = nPeaks, shape = peak.mean.shape,
                scale = peak.mean.scale
            )
        },
        lngamma = {
            peak.mean.pi <- simPICget(object, "peak.mean.pi")
            peak.mean.shape <- simPICget(object, "peak.mean.shape")
            peak.mean.rate <- simPICget(object, "peak.mean.rate")
            peak.mean.meanlog <- simPICget(object, "peak.mean.meanlog")
            peak.mean.sdlog <- simPICget(object, "peak.mean.sdlog")
            peak.means <- rlngamma(
                n = nPeaks, pi = peak.mean.pi, shape = peak.mean.shape,
                rate = peak.mean.rate, meanlog = peak.mean.meanlog,
                sdlog = peak.mean.sdlog
            )
        },
        stop("Invalid distribution: ", pm.distr)
    )

    peak.means.normalised <- peak.means / sum(peak.means)
    rowData(sim)$exp.peakmean <- peak.means.normalised
    return(sim)
}

#' Simulate true counts.
#'
#' Counts are simulated from a poisson distribution where each peak has a
#' mean, expected library size and proportion of accessible chromatin.
#'
#' @param sim SingleCellExperiment object containing simulation parameters.
#' @param object simPICcount object with simulation parameters.
#' @return SingleCellExperiment object with simulated true counts.
#' @importFrom SummarizedExperiment rowData colData<- colData colData<-
#' @importFrom SummarizedExperiment assays assays<-
#' @importFrom stats rbinom rpois

simPICsimulateTrueCounts <- function(object, sim) {
    nCells <- simPICget(object, "nCells")
    nPeaks <- simPICget(object, "nPeaks")
    peak.mean <- rowData(sim)$exp.peakmean
    lib.size <- colData(sim)$exp.libsize
    sparsity <- simPICget(object, "sparsity")

    if (nPeaks > length(sparsity)) {
        sparsity.max <- max(sparsity)
        sparsity.min <- min(sparsity)
        sparsity <- sample(seq(sparsity.min, sparsity.max), nPeaks, 
                                replace = TRUE)
    }

    true_counts <- matrix(nrow = nPeaks, ncol = nCells)

    sim_true_counts <- function(nCells) {
        lambda <- peak.mean * lib.size[nCells]
        count_vec <- rpois(as.numeric(nPeaks), lambda)
        zero_prop <- sparsity[nPeaks]

        if (length(count_vec[count_vec == 0]) == 0) {
            num_zero <- rep(0, nPeaks)
        } else {
            if (sum(count_vec == 0) > 0) {
                num_zero <- rbinom(
                    nPeaks, 1,
                    zero_prop - mean(count_vec[count_vec > 0] == 0)
                )
            } else {
                num_zero <- rbinom(nPeaks, 1, zero_prop)
            }
        }
        count_vec <- count_vec * num_zero
        count_vec
    }

    true_counts <- as(
                    sapply(seq_len(nCells),
                    sim_true_counts),
                    "dgCMatrix"
    )
    assays(sim, withDimnames = FALSE)$counts <- true_counts
    return(sim)
}
