#' Estimate simPIC simulation parameters
#'
#' Estimate simulation parameters for library size, peak means, and sparsity
#' from a real peak-by-cell input matrix.
#'
#' @param counts either a sparse peak by cell count matrix, or a
#'        SingleCellExperiment object containing count data to estimate
#'        parameters.
#' @param object simPICcount object to store estimated parameters and
#'        counts.
#' @param pm.distr statistical distribution for estimating peak mean
#'        parameters. Available distributions: gamma, weibull, lngamma, pareto.
#'        Default is lngamma.
#' @param method Simulation mode. Use \code{"single"} to estimate parameters
#'        for one cell type or \code{"groups"} for distinct cell types.
#' @param verbose logical variable. Prints the simulation progress if TRUE.
#'
#' @return simPICcount object containing all estimated parameters.
#' @examples
#' counts <- readRDS(system.file("extdata", "test.rds", package = "simPIC"))
#' est <- newsimPICcount()
#' est <- simPICestimate(counts, pm.distr = "lngamma")
#' @export
simPICestimate <- function(counts,
                        object = newsimPICcount(),
                        pm.distr = c(
                        "lngamma", "gamma", "weibull", "pareto"
                        ),
                        method = c("single","groups"),
                        verbose = TRUE) {
    UseMethod("simPICestimate")
}

#' @rdname simPICestimate
#' @importFrom methods as
#' @export
simPICestimate.SingleCellExperiment <- function(counts,
                                                object = newsimPICcount(),
                                                pm.distr = "lngamma",
                                                method = "single",
                                                verbose = TRUE) {
    checkmate::assert_choice(pm.distr, c(
        "gamma", "weibull",
        "pareto", "lngamma"
    ))
    checkmate::assert_choice(method, c("single", "groups"))
    counts <- getCounts(counts)
    simPICestimate(
        counts,
        object = object,
        pm.distr = pm.distr,
        method = method,
        verbose = verbose
    )
}

#' @rdname simPICestimate
#' @importFrom stats median
#' @importFrom Matrix colSums
#' @export
simPICestimate.dgCMatrix <- function(counts,
                            object = newsimPICcount(),
                            pm.distr = "lngamma",
                            method = "single",
                            verbose = TRUE) {
    checkmate::assertClass(object, "simPICcount")
    checkmate::assert_choice(pm.distr, c(
        "gamma", "weibull",
        "pareto", "lngamma"
    ))
    checkmate::assert_choice(method, c("single", "groups"))
    
    object <- setsimPICparameters(object,
        nPeaks = nrow(counts),
        nCells = ncol(counts),
        batchCells = ncol(counts),
        pm.distr = pm.distr
    )

    counts <- counts[, which(colSums(counts) != 0), drop=FALSE]

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
    object <- simPICestimatePeakMean(norm.counts, object, pm.distr, verbose)
    
    if(method == "groups")
    {
        message("estimating BCV...")
        object <- simPICEstBCV(counts, object, verbose)
    }
    return(object)
}

#' Estimate simPIC library size parameters
#'
#' Estimate the library size parameters for simPIC simulation.
#'
#' @param counts count matrix.
#' @param object simPICcount object to store estimated values.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @details
#' Parameters for the lognormal distribution are estimated by fitting the
#' library sizes using \code{\link[fitdistrplus]{fitdist}}. All the fitting
#' methods are tried and the fit with the best Cramer-von Mises statistic is
#' selected.
#'
#' @return simPICcount object with estimated library size parameters.
#'
#' @importFrom Matrix colSums
simPICestimateLibSize <- function(counts, object, verbose) {
    lib.size <- colSums(counts)

    fit <- selectFit(lib.size, "lnorm", verbose = verbose)

    lib.size.meanlog <- unname(fit$estimate["meanlog"])
    lib.size.sdlog <- unname(fit$estimate["sdlog"])

    object <- setsimPICparameters(object,
        lib.size.meanlog = lib.size.meanlog,
        lib.size.sdlog = lib.size.sdlog
    )
    return(object)
}

#' Select fit
#'
#' Trying two fitting methods and selecting the best one.
#'
#' @param data The data to fit.
#' @param distr Name of the distribution to fit.
#' @param verbose Logical. Whether to print progress messages.
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
    }
    return(fits[[selected]])
}

#' Estimate sparsity
#'
#' This function estimates cell sparsity from a normalized count matrix and
#' updates the parameters of a simPIC object accordingly.
#'
#'
#' @param norm.counts A normalized count matrix to estimate parameters from.
#' @param object simPICcount object to store estimated parameters.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @return simPICcount object with updated sparsity parameter.
#'
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
#' @param norm.counts Library-size normalized count matrix.
#' @param object simPICcount object to store estimated values.
#' @param pm.distr distribution parameter for peak means.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @details
#' Parameters for gamma distribution are estimated by fitting the mean
#' normalized counts using \code{\link[fitdistrplus]{fitdist}}.
#' All the fitting methods are tried and the fit with the best Cramer-von
#' Mises statistic is selected.
#' @return simPICcount object containing all estimated parameters
#' @importFrom Matrix rowMeans
#' @importFrom stats sd
simPICestimatePeakMean <- function(norm.counts, object, pm.distr, verbose) {
    logical_matrix <- norm.counts != 0
    norm.counts <- norm.counts * logical_matrix
    peak.means <- rowMeans(norm.counts)
    peak.means <- peak.means[peak.means != 0]
    message("using ", pm.distr, " distribution for estimating peak mean")
    switch(pm.distr,
        gamma = {
            fit <- selectFit(peak.means, "gamma", verbose = verbose)
            peak.mean.shape <- unname(fit$estimate["shape"])
            peak.mean.rate <- unname(fit$estimate["rate"])
            object <- setsimPICparameters(object,
                peak.mean.shape = peak.mean.shape,
                peak.mean.rate = peak.mean.rate
            )
        },
        weibull = {
            pseudomeans <- peak.means / sd(peak.means)
            fit <- selectFit(pseudomeans, "weibull", verbose = verbose)
            peak.mean.shape <- unname(fit$estimate["shape"])
            peak.mean.scale <- unname(fit$estimate["scale"] * sd(peak.means))
            object <- setsimPICparameters(object,
                peak.mean.shape = peak.mean.shape,
                peak.mean.scale = peak.mean.scale
            )
        },
        pareto = {
            fit <- with_fitdist_bindings(
                list(
                    dpareto = function(x, shape, scale, log = FALSE) {
                        actuar::dpareto(
                            x,
                            shape = shape,
                            scale = scale,
                            log = log
                        )
                    },
                    ppareto = function(q, shape, scale, lower.tail = TRUE,
                                       log.p = FALSE) {
                        actuar::ppareto(
                            q,
                            shape = shape,
                            scale = scale,
                            lower.tail = lower.tail,
                            log.p = log.p
                        )
                    }
                ),
                fitdistrplus::fitdist(
                    peak.means,
                    "pareto",
                    start = list(shape = 1.3, scale = 0.05),
                    control = list(maxit = 1000)
                )
            )
            peak.mean.shape <- unname(fit$estimate["shape"])
            peak.mean.scale <- unname(fit$estimate["scale"])
            object <- setsimPICparameters(object,
                peak.mean.shape = peak.mean.shape,
                peak.mean.scale = peak.mean.scale
            )
        },
        lngamma = {
            fit <- with_fitdist_bindings(
                list(
                    dlngamma = function(x, pi, shape, rate, meanlog, sdlog) {
                        get("dlngamma", envir = environment(simPICestimatePeakMean))(
                            x,
                            pi = pi,
                            shape = shape,
                            rate = rate,
                            meanlog = meanlog,
                            sdlog = sdlog
                        )
                    },
                    plngamma = function(q, pi, shape, rate, meanlog, sdlog) {
                        get("plngamma", envir = environment(simPICestimatePeakMean))(
                            q,
                            pi = pi,
                            shape = shape,
                            rate = rate,
                            meanlog = meanlog,
                            sdlog = sdlog
                        )
                    }
                ),
                fitdistrplus::fitdist(
                    peak.means,
                    "lngamma",
                    optim.method = "BFGS",
                    start = list(
                        pi = 0, shape = 0, rate = 0, meanlog = 0, sdlog = 1
                    ),
                    control = list(maxit = 1000)
                )
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
        },
        stop("Invalid distribution: ", pm.distr)
    )
    return(object)
}

with_fitdist_bindings <- function(bindings, expr) {
    eval_expr <- substitute(expr)
    search_name <- paste0(
        "simPIC_fitdist_",
        as.integer(stats::runif(1, min = 1, max = 1e9))
    )
    base_attach <- get("attach", envir = baseenv())
    base_attach(
        list2env(bindings, parent = emptyenv()),
        name = search_name,
        warn.conflicts = FALSE
    )
    on.exit({
        search_pos <- match(search_name, search())
        if (!is.na(search_pos)) {
            detach(pos = search_pos)
        }
    }, add = TRUE)
    eval(eval_expr, envir = parent.frame())
}

#' Estimate simPIC Biological Coefficient of Variation parameters
#'
#' Parameters are estimated using the \code{\link[edgeR]{estimateDisp}} function
#' in the \code{edgeR} package.
#'
#' @param counts counts matrix to estimate parameters from.
#' @param object simPICcount object to store estimated values in.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @details
#' The \code{\link[edgeR]{estimateDisp}} function is used to estimate the common
#' dispersion and prior degrees of freedom. See
#' \code{\link[edgeR]{estimateDisp}} for details. When estimating parameters on
#' simulated data we found a broadly linear relationship between the true
#' underlying common dispersion and the \code{edgeR} estimate, therefore we
#' apply a small correction, \code{disp = -0.3 + 0.15 * edgeR.disp}.
#'
#' @return simPICcount object with estimated values.

simPICEstBCV <- function(counts, object, verbose) {
    # Add dummy design matrix to avoid print statement
    design <- matrix(1, ncol(counts), 1)
    disps <- edgeR::estimateDisp(counts, design = design)
    bcv.common <- -0.3 + 0.15 * disps$common.dispersion
    
    bcv.df <- disps$prior.df
    
    object <- setsimPICparameters(
        object,
        bcv.common = bcv.common,
        bcv.df = bcv.df 
    )
    
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
