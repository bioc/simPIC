    #' Estimate simPIC simulation parameters
    #'
    #' Estimate parameters library size, peak means, and sparsity
    #' for simPIC simulation from a real peak by cell input matrix
    #'
    #' @param counts either a sparse peak by cell count matrix, or a
    #' SingleCellExperiment object containing a count matrix to estimate
    #' parameters.
    #' @param object simPICcount object to store estimated parameters and
    #' counts.
    #' @param pm.distr statistical distribution for estimating peak mean
    #' parameters.
    #' @param verbose logical variable. Prints the simulation progress if TRUE.
    #'
    #' @return simPICcount object containing all estimated parameters.
    #' @examples
    #' counts <- readRDS(system.file("extdata","test.rds", package = "simPIC"))
    #' est <- newsimPICcount()
    #' est <- simPICestimate(counts, pm.distr = "weibull")
    #'
    #' @export
    simPICestimate <- function(counts,
                        object = newsimPICcount(),
                        pm.distr = c(
                        "gamma", "weibull", "pareto","lngamma"
                        ),
                        verbose = TRUE) {
    UseMethod("simPICestimate")
    }

    #' @rdname simPICestimate
    #' @importFrom methods as
    #' @export
    simPICestimate.SingleCellExperiment <- function(counts,
                                                    object = newsimPICcount(),
                                                    pm.distr = c(
                                                        "gamma", "weibull",
                                                        "pareto", "lngamma"
                                                    ),
                                                    verbose = TRUE) {
        counts <- getCounts(counts)
        simPICestimate.dgCMatrix(
            as(counts, "dgCMatrix"), object, pm.distr,
            verbose
        )
    }

    #' @rdname simPICestimate
    #' @importFrom stats median
    #' @importFrom Matrix colSums
    #' @export
    simPICestimate.dgCMatrix <- function(counts,
                                    object = newsimPICcount(),
                                    pm.distr = "weibull",
                                    verbose = TRUE) {
    checkmate::assertClass(object, "simPICcount")
    checkmate::assert_choice(pm.distr, c(
    "gamma", "weibull",
    "pareto", "lngamma"
        )
    )

        object <- setsimPICparameters(object,
            nPeaks = nrow(counts),
            nCells = ncol(counts)
        )

        counts <- counts[, which(colSums(counts) != 0)]

        lib.sizes <- colSums(counts)
        lib.med <- median(lib.sizes)
        norm.counts <- as(
            t(t(as.matrix(counts)) / lib.sizes * lib.med),
            "dgCMatrix"
        )

        if (dim(counts)[1] > 0) {
            object <- setsimPICparameters(object)
        }

        if (verbose) {
            message("simPIC is:")
        }

        if (verbose) {
            message("estimating library size parameters...")
        }
        object <- simPICestimateLibSize(counts, object, verbose)

        if (verbose) {
            message("estimating sparsity...")
        }
        object <- simPICestimateSparsity(counts, object, verbose)

        if (verbose) {
            message("estimating peak mean parameters...")
        }
        object <-
            simPICestimatePeakMean(norm.counts, object, pm.distr, verbose)

        return(object)
    }

    #' Estimate simPIC library size parameters.
    #'
    #' Estimate the library size parameters for simPIC simulation.
    #'
    #' @param counts count matrix.
    #' @param object simPICcount object to store estimated values.
    #' @param verbose logical. To print messages or not.
    #'
    #' @details
    #' Parameters for the lognormal distribution are estimated by fitting the
    #' library sizes using \code{\link[fitdistrplus]{fitdist}}. All the fitting
    #' methods are tried and the fit with the best Cramer-von Mises statistic is
    #' selected.
    #'
    #' @return simPICcount object with estimated library size parameters
    #'
    #' @importFrom Matrix colSums
    simPICestimateLibSize <- function(counts, object, verbose) {
        lib.size <- colSums(counts)

        fit <- selectFit(lib.size, "lnorm", verbose = verbose)

        lib.size.meanlog <- unname(fit$estimate["meanlog"])
        lib.size.sdlog <- unname(fit$estimate["sdlog"])

        object <- setsimPICparameters(object,
            lib.size.meanlog = lib.size.meanlog,
            lib.size.sdlog = lib.size.sdlog,
            #lib.sizes = lib.size
        )
        return(object)
    }

    #' Select fit
    #'
    #' Trying a variety of fitting methods and selecting the best one
    #'
    #' @param data The data to fit
    #' @param distr Name of the distribution to fit
    #' @param verbose logical. To print messages or not
    #'
    #' @details
    #' The distribution is fitted to the data using each of the
    #' \code{\link[fitdistrplus]{fitdist}} fitting methods. The fit with the
    #' smallest Cramer-von Mises statistic is selected.
    #'
    #' @return The selected fit object

    selectFit <- function(data, distr, verbose = TRUE) {
        checkmate::assertNumeric(data, finite = TRUE, any.missing = FALSE)
        checkmate::assertString(distr)
        checkmate::assertFlag(verbose)

        # Sink output that sometimes happens when fitting
        sink(tempfile())
        on.exit(sink())

        fits <- list()

        try(
            fits$`MLE` <- fitdistrplus::fitdist(data, distr, method = "mle"),
            silent = TRUE
        )

        try(
            fits$`MGE (CvM)` <- fitdistrplus::fitdist(data, distr,
                method = "mge",
                gof = "CvM"
            ),
            silent = TRUE
        )

        scores <- fitdistrplus::gofstat(fits)$cvm

        # Flatten in case scores is a list, selecting the score with min cvm
        scores.flat <- unlist(scores)
        selected <- which(scores.flat == min(scores.flat, na.rm = TRUE))

        if (verbose) {
            # Work around to get name in case scores is a list
            name <- names(fits)[names(scores) == names(scores.flat)[selected]]
            # message("Selected ", name, " fit")
        }
        return(fits[[selected]])
    }

    #' Estimate simPIC peak sparsity.
    #'
    #' Extract the accessibility proportion (sparsity) of each cell among all
    #' peaksvfrom the input count matrix.
    #'
    #' @param norm.counts A sparse count matrix to estimate parameters from.
    #' @param object simPICcount object to store estimated parameters.
    #' @param verbose logical. To print messages or not.
    #'
    #' @return simPICcount object with updated non-zero cell proportion
    #' parameter.
    #'
    #' @details
    #' Vector of non-zero cell proportions of peaks is calculated by
    #' dividing the number of non-zero entries over the number of all cells
    #' for each peak.
    #' @importFrom Matrix rowSums

    simPICestimateSparsity <- function(norm.counts, object, verbose) {
        sparsity <- ifelse(rowMeans(norm.counts == 0) < 1,
            rowMeans(norm.counts == 0), 0
        )
        object <- setsimPICparameters(object, sparsity = sparsity)
        return(object)
    }

    #' Estimate simPIC peak means
    #'
    #' Estimate peak mean parameters for simPIC simulation
    #'
    #' @param norm.counts library size normalised counts matrix.
    #' @param object simPICcount object to store estimated values.
    #' @param pm.distr distribution parameter for peak means
    #' @param verbose logical. To print progress messages or not.
    #'
    #' @details
    #' Parameters for gamma distribution are estimated by fitting the mean
    #' normalised counts using \code{\link[fitdistrplus]{fitdist}}.
    #' All the fitting methods are tried and the fit with the best Cramer-von
    #' Mises statistic is selected.
    #' @return simPICcount object containing all estimated parameters
    #' @importFrom Matrix rowMeans
    #' @importFrom stats sd
    #' @importFrom actuar dpareto ppareto
    simPICestimatePeakMean <- function(norm.counts, object, pm.distr, verbose) {
        #lib.sizes <- simPICget(object, "lib.sizes")
        logical_matrix <- norm.counts != 0
        norm.counts <- norm.counts * logical_matrix
        peak.means <- rowMeans(norm.counts)
        peak.means <- peak.means[peak.means != 0]
        if (pm.distr == "gamma") {
            message("using gamma distribution for estimating peak mean")
            fit <- selectFit(peak.means, "gamma", verbose = verbose)
            peak.mean.shape <- unname(fit$estimate["shape"])
            peak.mean.rate <- unname(fit$estimate["rate"])
            object <- setsimPICparameters(object,
                peak.mean.shape = peak.mean.shape,
                peak.mean.rate = peak.mean.rate
            )
        } else if (pm.distr == "weibull") {
            message("using weibull distribution for estimating peak mean")
            pseudomeans <- peak.means / sd(peak.means)
            fit <- selectFit(pseudomeans, "weibull", verbose = verbose)
            peak.mean.shape <- unname(fit$estimate["shape"])
            peak.mean.scale <- unname(fit$estimate["scale"] * sd(peak.means))
            object <- setsimPICparameters(object,
                peak.mean.shape = peak.mean.shape,
                peak.mean.scale = peak.mean.scale
            )
        } else if (pm.distr == "pareto") {
            message("using pareto distribution for estimating peak mean")
            fit <- fitdistrplus::fitdist(peak.means, "pareto",
                start = list(shape = 1.3, scale = 0.05), control = list(maxit = 1000)
            )
            peak.mean.shape <- unname(fit$estimate["shape"])
            peak.mean.scale <- unname(fit$estimate["scale"])
            object <- setsimPICparameters(object,
                peak.mean.shape = peak.mean.shape,
                peak.mean.scale = peak.mean.scale
            )
        } else if (pm.distr == "lngamma") {
            message("using lngamma distribution for estimating peak mean")
            fit <- fitdistrplus::fitdist(peak.means, "lngamma",
                optim.method = "BFGS",
                start = list(pi = 0, shape = 0, rate = 0, meanlog = 0, sdlog = 1),
                control = list(maxit = 1000)
            )
            peak.mean.pi <- unname(fit$estimate["pi"])
            peak.mean.shape <- unname(fit$estimate["shape"])
            peak.mean.rate <- unname(fit$estimate["rate"])
            peak.mean.meanlog <- unname(fit$estimate["meanlog"])
            peak.mean.sdlog <- unname(fit$estimate["sdlog"])
            object <- setsimPICparameters(object,
                peak.mean.pi = peak.mean.pi,
                peak.mean.shape = peak.mean.shape, 
                peak.mean.rate = peak.mean.rate,
                peak.mean.meanlog = peak.mean.meanlog,
                peak.mean.sdlog = peak.mean.sdlog
            )
        } else {
            stop("Invalid distribution choice")
        }
        return(object)
    }
    dlngamma <- function(x, pi, shape, rate, meanlog, sdlog) {
        pi <- zero_one(pi)
        shape <- positive(shape)
        rate <- positive(rate)
        sdlog <- abs(sdlog)
        pi * dgamma(x, shape, rate) + (1 - pi) * dlnorm(x, meanlog, sdlog)
    }
    plngamma <- function(q, pi, shape, rate, meanlog, sdlog) {
        pi <- zero_one(pi)
        shape <- positive(shape)
        rate <- positive(rate)
        sdlog <- abs(sdlog)
        pi * pgamma(q, shape, rate) + (1 - pi) * plnorm(q, meanlog, sdlog)
    }
    rlngamma <- function(n, pi, shape, rate, meanlog, sdlog) {
        pi <- zero_one(pi)
        shape <- positive(shape)
        rate <- positive(rate)
        sdlog <- abs(sdlog)
        s1 <- rgamma(n, shape, rate)
        s2 <- rlnorm(n, meanlog, sdlog)
        ind <- runif(n) > pi
        s1[ind] <- s2[ind]
        s1
    }
    zero_one <- \(x) 1 / (1 + exp(-x))
    positive <- \(x) exp(x)
