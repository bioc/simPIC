#' simPIC simulation
#'
#' Simulate a peak-by-cell count matrix using simPIC methods.
#'
#' @param object simPICcount object with simulation parameters.
#' See \code{\link{simPICcount}} for details.
#' @param pm.distr distribution parameter for peak means.
#' Available distributions: gamma, weibull, lngamma, pareto.
#' Default is lngamma.
#' @param method Simulation mode. Use \code{"single"} to simulate one cell type
#'        or \code{"groups"} to simulate distinct cell types.
#' @param verbose Logical. Whether to print progress messages.
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
#'  The simulation involves the following steps:
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
#' \code{\link[SummarizedExperiment]{colData}} (for cell-specific information),
#' \code{\link[SummarizedExperiment]{rowData}} (for peak-specific information),
#' or \code{\link[SummarizedExperiment]{assays}} (for peak-by-cell matrices).
#' @examples
#' # default simulation
#' sim <- simPICsimulate(pm.distr = "lngamma")
#' @importFrom SummarizedExperiment rowData rowData<-
#' @importFrom SummarizedExperiment colData colData<-
#' @importFrom SummarizedExperiment assays assays<-
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @export
simPICsimulate <- function(object = newsimPICcount(),
                    pm.distr = "lngamma",
                    method = c("single","groups"),
                    verbose = TRUE,
                    ...) {
    checkmate::assertClass(object, "simPICcount")
    
    method <- match.arg(method)
    
    if (verbose) {
        message("simPIC is:")
        message("updating parameters...")
    }
    object <- setsimPICparameters(object, pm.distr = pm.distr, ...)
    validObject(object)
    # Set random seed
    seed <- simPICget(object, "seed")
    # Get the required parameters
    nCells <- simPICget(object, "nCells")
    nPeaks <- simPICget(object, "nPeaks")
    nGroups <- simPICget(object, "nGroups")
    group.prob <- simPICget(object, "group.prob")
    nBatches <- simPICget(object, "nBatches")
    batch.cells <- simPICget(object, "batchCells")
    
    # Ensure batch.cells matches nCells
    if (sum(batch.cells) != nCells) {
        warning("Sum of batchCells (", sum(batch.cells), 
                ") does not equal nCells (", nCells, "). Adjusting automatically.")
        batch.cells <- rep(floor(nCells / nBatches), nBatches)
        remainder <- nCells - sum(batch.cells)
        if (remainder > 0) {
            batch.cells[seq_len(remainder)] <- batch.cells[seq_len(remainder)] + 1
        }
        object <- setsimPICparameters(object, "batchCells" = batch.cells)
        object <- setsimPICparameters(object, "nCells" = sum(batch.cells))
        nCells <- sum(batch.cells)
    }
    
    if (nGroups == 1 && method == "groups") {
        warning("nGroups is 1, switching to default mode")
        method <- "single"
    }
    
    if (verbose) {
        message("setting up SingleCellExperiment object...")
    }
    
    cell.names <- paste0("Cell", seq_len(nCells))
    peak.names <- paste0("Peak", seq_len(nPeaks))
    batch.names <- paste0("Batch", seq_len(nBatches))
    if (method == "groups") {
        group.names <- paste0("Group", seq_len(nGroups))
    } else {
        paste0("no method")
    }
    # Create SingleCellExperiment to store simulation
    cells <- data.frame(Cell = cell.names)
    rownames(cells) <- cell.names
    peaks <- data.frame(Peak = peak.names)
    rownames(peaks) <- peak.names

    sim <- SingleCellExperiment(
        rowData = peaks, colData = cells, metadata = list(Params = object)
    )
    
    object <- setsimPICparameters(object, ...)
    validObject(object)
    # Make batches vector which is the index of object$batch.cells repeated
    # object$batch.cells[index] times
    batches <- lapply(seq_len(nBatches), function(i, b) {
        rep(i, b[i])
    }, b =  batch.cells)
    batches <- unlist(batches)
    colData(sim)$Batch <- batch.names[batches]
    
    withr::with_seed(seed, {
        if (method != "single") {
            groups <- sample(seq_len(nGroups), nCells,
                    prob = group.prob,
                    replace = TRUE
            )
            colData(sim)$Group <- factor(
                group.names[groups],
                levels = group.names
            )
        }
    if (verbose) {
        message("Simulating library size...")
    }
    sim <- simPICsimulateLibSize(object, sim)
        gc()
    if (verbose) {
        message("Simulating peak mean...")
    }
    sim <- simPICsimulatePeakMean(object, sim, pm.distr, verbose)
    
    if (nBatches > 1) {
            if (verbose) {
                message("Simulating batch effects...")
            }
            sim <- simPICsimBatchEffects(object, sim)
        }
        gc()
        sim <- simPICsimBatchCellMeans(object, sim)
        if (method == "single") {
            sim <- simPICsimSingleCellMeans(object, sim)
        } else if (method == "groups"){
            message("Simulating group DA")
            sim <- simPICsimulatemultiDA(object, sim)
            gc()
            message("Simulating group (cell) means")
            sim <- simPICsimulateGroupCellMeans(object, sim)
            message("Simulating BCV")
            sim <- simPICsimulateBCVMeans(object,sim)
        }
        
        if(method == "single")
        {
            message("Simulating true counts...")
            sim <- simPICsimulateTrueCounts(object, sim)
            rownames(BiocGenerics::counts(sim)) <- peak.names
            colnames(BiocGenerics::counts(sim)) <- cell.names
        } else if(method == "groups")
        {
            message("Simulating true counts groups...")
            sim <- simPICsimulateTrueCountsGroups(object,sim)
            rownames(BiocGenerics::counts(sim)) <- peak.names
            colnames(BiocGenerics::counts(sim)) <- cell.names
        }
        else 
        {
            message("invalid method, provide single or groups")
        }
    })
    if (verbose) {
        message("Done!!")
    }
    
    return(sim)
}

#' @rdname simPICsimulate
#' @export
simPICsimulatesingle <- function(object = newsimPICcount(),
                                verbose = TRUE, ...){
    sim <- simPICsimulate(
        object = object,
        method = "single",
        verbose = verbose,
        ...
    )
    return(sim)
}


#' @rdname simPICsimulate
#' @export
simPICsimulatemulti <- function(object = newsimPICcount(),
                                pm.distr = "lngamma", 
                                method = c("groups"),
                                verbose = TRUE, ...) {
    sim <- simPICsimulate(
        object = object,
        pm.distr = pm.distr,
        method = "groups",
        verbose = verbose,
        ...
    )
    return(sim)
}


#' Simulate simPIC library sizes
#'
#' Generate library sizes for cells in simPIC simulation based on the
#' estimated values of mus and sigmas.
#'
#' @param sim SingleCellExperiment object containing simulation parameters.
#' @param object simPICcount object with simulation parameters.
#' @param verbose Logical. Whether to print progress messages.
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

#' Simulate simPIC peak means
#'
#' Generate peak means for cells in simPIC simulation based on the estimated
#' values of shape and rate parameters.
#'
#' @param sim SingleCellExperiment object containing simulation parameters.
#' @param object simPICcount object with simulation parameters.
#' @param pm.distr distribution parameter for peak means.
#' Available distributions: gamma, weibull, lngamma, pareto.
#' Default is lngamma.
#' @param verbose logical. Whether to print progress messages.
#'
#' @return SingleCellExperiment object with simulated peak means.
#'
#' @importFrom SummarizedExperiment rowData rowData<-
#' @importFrom stats rgamma rweibull dgamma dlnorm pgamma plnorm runif
#' @importFrom Matrix sparseMatrix
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
            peak.means <- actuar::rpareto(
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

#' Simulate batch effects
#'
#' Simulate batch effects. Batch effect factors for each batch are produced
#' using \code{\link{getLNormFactors}} and these are added along with updated
#' means for each batch.
#'
#' @param sim SingleCellExperiment to add batch effects to.
#' @param object simPICcount object with simulation parameters.
#'
#' @return SingleCellExperiment with simulated batch effects.
#'
#' @importFrom SummarizedExperiment rowData rowData<-
#'
#' @keywords internal
simPICsimBatchEffects <- function(object,sim) {
    nPeaks <- simPICget(object, "nPeaks")
    nBatches <- simPICget(object, "nBatches")
    batch.facLoc <- simPICget(object, "batch.facLoc")
    batch.facScale <- simPICget(object, "batch.facScale")
    batch.rmEffect <- simPICget(object, "batch.rmEffect")
    peak.means.normalised <- rowData(sim)$exp.peakmean
    
    batch.facLoc <- rep(batch.facLoc, nBatches)
    batch.facScale <- rep(batch.facScale, nBatches)
    
    for (idx in seq_len(nBatches)) {
        batch.facs <- getLNormFactors(
            nPeaks, 1, 0.5, batch.facLoc[idx], batch.facScale[idx]
        )
        
        if (batch.rmEffect) {
            batch.facs <- rep(1, length(batch.facs))
        }
        
        rowData(sim)[[paste0("BatchFacBatch", idx)]] <- batch.facs
    }
    
    return(sim)
}

#' Simulate batch means
#'
#' Simulate a mean for each peak in each cell incorporating batch effect
#' factors.
#'
#' @param sim SingleCellExperiment to add batch means to.
#' @param object simPICcount object with simulation parameters.
#'
#' @return SingleCellExperiment with simulated batch means.
#'
#' @importFrom SummarizedExperiment rowData rowData<- colData
#'
#' @keywords internal
simPICsimBatchCellMeans <- function(object, sim) {
    nBatches <- simPICget(object, "nBatches")
    cell.names <- colData(sim)$Cell
    peak.names <- rowData(sim)$Peak
    peak.means <- rowData(sim)$exp.peakmean
    
    if (nBatches > 1) {
        batches <- colData(sim)$Batch
        batch.names <- unique(batches)
        
        batch.facs.peak <- as.matrix(
            rowData(sim)[, paste0("BatchFac", batch.names)]
        )
        batch.facs.cell <- as.matrix(
            batch.facs.peak[, as.numeric(factor(batches))]
        )
    } else {
        nCells <- length(colData(sim)$Cell)
        nPeaks <- simPICget(object, "nPeaks")
        
        batch.facs.cell <- matrix(1, ncol = nCells, nrow = nPeaks)
    }
    
    batch.means.cell <- batch.facs.cell * peak.means
    
    colnames(batch.means.cell) <- cell.names
    rownames(batch.means.cell) <- peak.names
    assays(sim)$BatchCellMeans <- batch.means.cell
    
    return(sim)
}


#' Simulate group differential accessibility
#'
#' Simulate differential accessibility. Differential accessibility factors for each
#' group are produced using \code{\link{getLNormFactors}} and these are added
#' along with updated means for each group. For paths care is taken to make sure
#' they are simulated in the correct order.
#'
#' @param sim SingleCellExperiment to add differential accessibility to.
#' @param object simPICcount object with simulation parameters.
#'
#' @return SingleCellExperiment with simulated differential accessibility.
#'
#' @name simPICsimulatemultiDA
#' @rdname simPICsimulatemultiDA
#' @importFrom SummarizedExperiment rowData
simPICsimulatemultiDA <- function(object, sim) {
    nPeaks <- simPICget(object, "nPeaks")
    nGroups <- simPICget(object, "nGroups")
    da.prob <- simPICget(object, "da.prob")
    da.downProb <- simPICget(object, "da.downProb")
    da.facLoc <- simPICget(object, "da.facLoc")
    da.facScale <- simPICget(object, "da.facScale")
    peak.mean <- rowData(sim)$exp.peakmean
    
    ## update
    da.prob <- rep(da.prob, nGroups)
    da.downProb <- rep(da.downProb, nGroups)
    da.facLoc <- rep(da.facLoc, nGroups)
    da.facScale <- rep(da.facScale, nGroups)
    
    for (idx in seq_len(nGroups)) {
        da.facs <- getLNormFactors(
            nPeaks,
            da.prob[idx],
            da.downProb[idx],
            da.facLoc[idx],
            da.facScale[idx]
        )
        group.means.peak <- peak.mean * da.facs
        rowData(sim)[[paste0("DAFacGroup", idx)]] <- da.facs
    }
    return(sim)
}

#' Simulate cell means
#' 
#' Simulate a peak by cell matrix given the mean accessibility for each peak in
#' each cell. Cells start with the mean accessibility for the group they belong
#' to (when simulating groups). The selected means are adjusted for each cell's
#' expected library size.
#' 
#' @param sim SingleCellExperiment to add cell means to.
#' @param object simPIC object with simulation parameters.
#' 
#' @return SingleCellExperiment with added cell means.
#' 
#' @name simPICsimCellMeans
#' 
#' @rdname simPICsimCellMeans
#' @importFrom SummarizedExperiment rowData colData assays assays<-
simPICsimSingleCellMeans <- function(object,sim) {
    nCells <- length(colData(sim)$Cell)
    cell.names <- colData(sim)$Cell
    peak.names <- rowData(sim)$Peak
    exp.libsize <- colData(sim)$exp.libsize
    batch.means.cell <- assays(sim)$BatchCellMeans

    cell.means.peak <- batch.means.cell
    cell.props.peak <- t(t(cell.means.peak) / colSums(cell.means.peak))
    base.means.cell <- t(t(cell.props.peak) * exp.libsize)

    colnames(base.means.cell) <- cell.names
    rownames(base.means.cell) <- peak.names
    assays(sim)$BaseCellMeans <- base.means.cell

    return(sim)
}

#' @name simPICsimCellMeans
#' 
#' @rdname simPICsimCellMeans
#' @importFrom SummarizedExperiment rowData colData assays assays<-

simPICsimulateGroupCellMeans <- function(object, sim) {
    nGroups <-  simPICget(object, "nGroups")
    cell.names <- colData(sim)$Cell
    peak.names <- rowData(sim)$Peak
    groups <- colData(sim)$Group
    group.names <- levels(groups)
    exp.libsize <- colData(sim)$exp.libsize
    batch.means.cell <- assays(sim)$BatchCellMeans
    
    group.facs.peak <- rowData(sim)[, paste0("DAFac", group.names)]
    cell.facs.peak <- as.matrix(group.facs.peak[, paste0("DAFac", groups)])
    cell.means.peak <- batch.means.cell * cell.facs.peak 

    cell.props.peak <- t(t(cell.means.peak) / colSums(cell.means.peak))

    base.means.cell <- t(t(cell.props.peak) * exp.libsize) 
    
    colnames(base.means.cell) <- cell.names
    rownames(base.means.cell) <- peak.names
    assays(sim)$BaseCellMeans <- base.means.cell
    
    return(sim)
}

#' Get accessibility factors
#'
#' Randomly generate multiplication factors from a log-normal distribution.
#'
#' @param n.facs Number of factors to generate.
#' @param sel.prob Probability that a factor will be selected to be different
#'        from 1.
#' @param neg.prob Probability that a selected factor is less than one.
#' @param fac.loc Location parameter for the log-normal distribution.
#' @param fac.scale Scale factor for the log-normal distribution.
#'
#' @return Vector containing generated factors.
#'
#' @importFrom stats rbinom rlnorm
#'
#' @keywords internal
getLNormFactors <- function(n.facs, sel.prob, neg.prob, fac.loc, fac.scale) {
    # Validate inputs
    if (!is.numeric(n.facs) || n.facs <= 0 || round(n.facs) != n.facs) {
        stop("n.facs must be a positive integer")
    }
    if (!is.numeric(sel.prob) || sel.prob < 0 || sel.prob > 1) {
        stop("sel.prob must be a probability between 0 and 1")
    }
    if (!is.numeric(neg.prob) || neg.prob < 0 || neg.prob > 1) {
        stop("neg.prob must be a probability between 0 and 1")
    }
    if (!is.numeric(fac.loc)) {
        stop("fac.loc must be numeric")
    }
    if (!is.numeric(fac.scale)) {
        stop("fac.scale must be numeric")
    }
    
    is.selected <- as.logical(rbinom(n.facs, 1, sel.prob))
    n.selected <- sum(is.selected)
    dir.selected <- (-1)^rbinom(n.selected, 1, neg.prob)
    facs.selected <- rlnorm(n.selected, fac.loc, fac.scale)
    
    # Count and print the number of NAs created
    num_na <- sum(is.na(facs.selected))
    if (num_na > 0) {
        message(sprintf("Number of NAs generated: %d", num_na))
    }
    
    # Re-run rlnorm only on the NA values
    while (any(is.na(facs.selected))) {
        na.indices <- which(is.na(facs.selected))
        
        facs.selected[na.indices] <- rlnorm(length(na.indices), fac.loc, fac.scale)
    }
    # Reverse directions for factors that are less than one
    dir.selected[facs.selected < 1] <- -1 * dir.selected[facs.selected < 1]
    factors <- rep(1, n.facs)
    factors[is.selected] <- facs.selected^dir.selected
    
    return(factors)
}


#' Simulate BCV means
#'
#' Simulate means for each peak in each cell that are adjusted to follow a
#' mean-variance trend using Biological Coefficient of Variation taken from
#' and inverse gamma distribution.
#'
#' @param sim SingleCellExperiment to add BCV means to.
#' @param object simPICcount object with simulation parameters.
#'
#' @return SingleCellExperiment with simulated BCV means.
#'
#' @importFrom SummarizedExperiment rowData colData assays assays<-
#' @importFrom stats rchisq rgamma
simPICsimulateBCVMeans <- function(object,sim) {
    cell.names <- colData(sim)$Cell
    peak.names <- rowData(sim)$Peak
    nPeaks <- simPICget(object, "nPeaks")
    nCells <- length(colData(sim)$Cell)
    bcv.common <- simPICget(object, "bcv.common")
    bcv.df <- simPICget(object, "bcv.df")
    base.means.cell <- assays(sim)$BaseCellMeans
    if (is.finite(bcv.df)) {
        bcv <- (bcv.common + (1 / sqrt(base.means.cell))) *
            sqrt(bcv.df / rchisq(nPeaks, df = bcv.df))
    } else {
        warning("'bcv.df' is infinite. This parameter will be ignored.")
        bcv <- (bcv.common + (1 / sqrt(base.means.cell)))
    }
    means.cell <- matrix(
        rgamma(
            as.numeric(nPeaks) * as.numeric(nCells),
            shape = 1 / (bcv^2),
            scale = base.means.cell * (bcv^2)
        ),
        nrow = nPeaks, ncol = nCells
    )
    
    colnames(means.cell) <- cell.names
    rownames(means.cell) <- peak.names
    
    assays(sim)$BCV <- bcv
    assays(sim)$CellMeans <- means.cell
    
    return(sim)
}


#' Simulate true counts
#'
#' Counts are simulated from a Poisson distribution where each peak has a
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

        true_counts <- vapply(seq_len(nCells), 
                            sim_true_counts, 
                            integer(nPeaks))
    
    assays(sim, withDimnames = FALSE)$counts <- true_counts
    sim <- ensureCountsFirst(sim)
    return(sim)
}





#' Simulate true counts groups
#'
#' Counts are simulated from a Poisson distribution where each peak has a
#' mean, expected library size and proportion of accessible chromatin.
#'
#' @param sim SingleCellExperiment object containing simulation parameters.
#' @param object simPICcount object with simulation parameters.
#' @return SingleCellExperiment object with simulated true counts.
#' @importFrom SummarizedExperiment rowData colData<- colData colData<-
#' @importFrom SummarizedExperiment assays assays<-
#' @importFrom stats rbinom rpois

simPICsimulateTrueCountsGroups <- function(object, sim) {
    cell.names <- colData(sim)$Cell
    peak.names <- rowData(sim)$Peak
    nCells <- simPICget(object, "nCells")
    nPeaks <- simPICget(object, "nPeaks")
    sparsity <- simPICget(object, "sparsity")
    cell.means <- assays(sim)$CellMeans
    
    if (nPeaks > length(sparsity)) {
        sparsity.max <- max(sparsity)
        sparsity.min <- min(sparsity)
        sparsity <- sample(seq(sparsity.min, sparsity.max), nPeaks, 
                    replace = TRUE)
    }
    
    zero_prop <- sparsity[nPeaks]
    
    true.counts <- matrix(
        rpois(
            as.numeric(nPeaks)* as.numeric(nCells),
            lambda = cell.means
        ),
        nrow = nPeaks, ncol = nCells
    )
    
    num_zero <- rbinom(nPeaks, 1, zero_prop)
    
    true.counts <- true.counts * num_zero
    
    colnames(true.counts) <- cell.names
    rownames(true.counts) <- peak.names
    
    assays(sim,withDimnames = FALSE)$counts <- true.counts
    sim <- ensureCountsFirst(sim)
    return(sim)
    
}
