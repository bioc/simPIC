#' Run the Packaged Microglia splatPop Example
#'
#' Run a streamlined population-scale Microglia example using packaged
#' \code{SingleCellExperiment}, sample mapping, peak annotation, and VCF files.
#' The workflow follows the splatPop estimation pattern by estimating
#' population means from retained donor-batch units while estimating
#' single-cell behaviour from the donor-batch unit with the most cells.
#'
#' @param sce Optional \code{SingleCellExperiment} or path to an RDS file
#' containing the packaged Microglia SCE. If \code{NULL}, the packaged example
#' is used.
#' @param sample.map Optional data.frame or path to a tab-separated sample map
#' linking VCF donor IDs to Microglia donor IDs. If \code{NULL}, the packaged
#' example is used.
#' @param peak.annot Optional data.frame or path to the chr14 peak annotation
#' table. If \code{NULL}, the packaged example is used.
#' @param vcf Optional \code{CollapsedVCF} or path to the chr14 VCF used for
#' the Microglia example. If \code{NULL}, the packaged example is used.
#' @param params Optional \code{SplatPopParams}. If \code{NULL}, a new object is
#' created with \code{pop.cv.bins}.
#' @param min.cells Minimum number of cells required per donor-batch unit for
#' inclusion in the aggregated population means. Units must exceed this
#' threshold.
#' @param pop.cv.bins Number of CV bins used when \code{params} is created
#' internally.
#' @param eqtl.n Proportion of peaks to simulate with caQTL effects. The
#' packaged example uses a modest default to keep the workflow lightweight
#' while still including genetic effects.
#' @param similarity.scale Similarity scale passed to
#' \code{splatter::setParams()} before simulation.
#' @param eqtl.dist Maximum cis-distance used when assigning caQTL effects.
#' @param sparsify Logical. Whether to sparsify the simulated output.
#' @param pca.ntop Number of variable peaks used for PCA plots.
#' @param pca.components Number of principal components to compute for the PCA
#' summaries.
#' @param plot.n.samples Number of libraries to show in the default PCA
#' comparison panels. The top libraries are chosen by real-data cell count.
#' @param plot.samples Optional character vector of specific libraries to show
#' in the PCA comparison panels. When \code{NULL}, the top
#' \code{plot.n.samples} libraries are used.
#' @param comparison.batch Optional real-data batch/library label used for the
#' Poptrial-style comparison panels. When supplied, the real data are subset to
#' this batch and the simulated data are matched by the same sample IDs.
#' @param comparison.pathology Optional pathology label used together with
#' \code{comparison.batch} for Poptrial-style panels, for example
#' \code{"earlyAD"}.
#' @param seed Optional random seed.
#' @param verbose Logical. Whether to print progress messages.
#' @param ... Additional parameters passed to \code{splatter::setParams()} prior
#' to simulation.
#'
#' @return A list containing:
#' \describe{
#' \item{\code{sce}}{Filtered real-data Microglia \code{SingleCellExperiment}.}
#' \item{\code{bigcounts}}{The donor-batch subset with the most cells used for
#' count-based estimation.}
#' \item{\code{aggregated}}{Aggregated donor-batch means as a
#' \code{SingleCellExperiment}.}
#' \item{\code{params}}{Estimated \code{SplatPopParams} after Microglia-specific
#' parameter updates.}
#' \item{\code{sim}}{Simulated \code{SingleCellExperiment}.}
#' \item{\code{plot_samples}}{Libraries used for the default PCA comparison
#' panels.}
#' \item{\code{plots}}{A named list containing real and simulated PCA plot
#' summaries, side-by-side comparison panels, and bluster-based silhouette and
#' neighborhood-purity panels.}
#' \item{\code{sample_map}}{The cleaned sample map used for alignment.}
#' \item{\code{peak_annot}}{The chr14 peak annotation used for simulation.}
#' \item{\code{vcf}}{The aligned VCF used for simulation.}
#' }
#'
#' @examples
#' \donttest{
#' if (requireNamespace("splatter", quietly = TRUE) &&
#'     requireNamespace("VariantAnnotation", quietly = TRUE) &&
#'     requireNamespace("bluster", quietly = TRUE)) {
#' out <- simPICMicrogliaExample(verbose = FALSE)
#' out$plot_samples
#' names(out$plots$comparison)
#' }
#' }
#' @export
simPICMicrogliaExample <- function(sce = NULL,
                                   sample.map = NULL,
                                   peak.annot = NULL,
                                   vcf = NULL,
                                   params = NULL,
                                   min.cells = 20L,
                                   pop.cv.bins = 10L,
                                   eqtl.n = 0.05,
                                   similarity.scale = 2,
                                   eqtl.dist = 1e8,
                                   sparsify = FALSE,
                                   pca.ntop = 100L,
                                   pca.components = 5L,
                                   plot.n.samples = 4L,
                                   plot.samples = NULL,
                                   comparison.batch = NULL,
                                   comparison.pathology = NULL,
                                   seed = NULL,
                                   verbose = TRUE,
                                   ...) {
    simPICrequireNamespace("splatter", "simPICMicrogliaExample()")
    simPICrequireNamespace("VariantAnnotation", "simPICMicrogliaExample()")
    simPICrequireNamespace("scuttle", "simPICMicrogliaExample()")

    checkmate::assertFlag(verbose)
    checkmate::assertFlag(sparsify)
    checkmate::assertCount(min.cells, positive = TRUE)
    checkmate::assertCount(pop.cv.bins, positive = TRUE)
    checkmate::assertNumber(eqtl.n, lower = 0, upper = 1, finite = TRUE)
    checkmate::assertNumber(similarity.scale, lower = 0, finite = TRUE)
    checkmate::assertNumber(eqtl.dist, lower = 0, finite = TRUE)
    checkmate::assertCount(pca.ntop, positive = TRUE)
    checkmate::assertCount(pca.components, positive = TRUE)
    checkmate::assertCount(plot.n.samples, positive = TRUE)
    checkmate::assertCharacter(plot.samples, null.ok = TRUE, min.chars = 1)
    checkmate::assertString(comparison.batch, null.ok = TRUE, min.chars = 1)
    checkmate::assertString(comparison.pathology, null.ok = TRUE, min.chars = 1)
    checkmate::assertInt(seed, null.ok = TRUE)

    paths <- simPICmicrogliaExamplePaths()

    if (is.null(sce)) {
        sce <- readRDS(paths$sce)
    } else if (is.character(sce) && length(sce) == 1L) {
        sce <- readRDS(sce)
    }
    checkmate::assertClass(sce, "SingleCellExperiment")

    if (is.null(sample.map)) {
        sample.map <- paths$sample_map
    }
    if (is.null(peak.annot)) {
        peak.annot <- paths$peak_annot
    }
    if (is.null(vcf)) {
        vcf <- paths$vcf
    }

    vcf <- simPICreadMicrogliaVcf(vcf)

    if (!is.null(seed)) {
        withr::local_seed(seed)
    }

    if (verbose) {
        message("Preparing packaged Microglia example inputs...")
    }
    prepared <- simPICprepareMicrogliaExample(
        sce = sce,
        sample.map = sample.map,
        peak.annot = peak.annot,
        vcf = vcf,
        min.cells = min.cells,
        force.keep.batch = comparison.batch,
        force.keep.pathology = comparison.pathology
    )

    batch.cells <- as.integer(table(prepared$sce$sample_batch))
    if (ncol(prepared$vcf) != length(batch.cells)) {
        stop(
            "The aligned VCF sample count (", ncol(prepared$vcf),
            ") must match the retained donor-batch count (",
            length(batch.cells), ")."
        )
    }

    if (verbose) {
        message("Estimating Microglia population parameters...")
    }
    params <- tryCatch(
        splatPopEstimatePeak(
            counts = prepared$bigcounts,
            means = getCounts(prepared$aggregated),
            params = params,
            sample.col = "Sample",
            batch.col = "Batch",
            aggregate.by = "sample_batch",
            pop.cv.bins = pop.cv.bins,
            verbose = verbose
        ),
        error = function(e) {
            if (!grepl("negative 'bcv.common'", conditionMessage(e), fixed = TRUE)) {
                stop(e)
            }

            warning(
                "The direct simPIC BCV correction was invalid for the ",
                "Microglia bigcounts subset. Continuing with the estimated ",
                "splatPop parameters and pop.quant.norm = FALSE.",
                call. = FALSE
            )

            splatPopEstimatePeak(
                counts = prepared$bigcounts,
                means = getCounts(prepared$aggregated),
                params = params,
                sample.col = "Sample",
                batch.col = "Batch",
                aggregate.by = "sample_batch",
                pop.cv.bins = pop.cv.bins,
                apply.bcv.correction = FALSE,
                verbose = verbose
            )
        }
    )
    params <- splatter::setParams(params, pop.quant.norm = FALSE)
    params <- splatter::setParams(
        params,
        batchCells = batch.cells,
        batch.size = 1L,
        eqtl.n = eqtl.n,
        similarity.scale = similarity.scale,
        eqtl.dist = eqtl.dist,
        ...
    )

    if (verbose) {
        message("Simulating Microglia population data with splatPop...")
    }
    sim <- withCallingHandlers(
        splatPopSimulatePeak(
            params = params,
            vcf = prepared$vcf,
            gff = prepared$peak.annot,
            sparsify = sparsify,
            verbose = verbose
        ),
        warning = function(w) {
            if (grepl("No SNP within eqtl.dist limit", conditionMessage(w), fixed = TRUE)) {
                invokeRestart("muffleWarning")
            }
        }
    )

    sim.cd <- SummarizedExperiment::colData(sim)
    if (!("Sample" %in% colnames(sim.cd))) {
        sample.candidates <- c("Individual", "sample", "ind", "Donor")
        match.idx <- sample.candidates[sample.candidates %in% colnames(sim.cd)]
        if (length(match.idx) > 0) {
            SummarizedExperiment::colData(sim)$Sample <- sim.cd[[match.idx[1]]]
        }
    }
    if (!("Batch" %in% colnames(sim.cd))) {
        batch.candidates <- c("batch", "Batch", "Experiment")
        match.idx <- batch.candidates[batch.candidates %in% colnames(sim.cd)]
        if (length(match.idx) > 0) {
            SummarizedExperiment::colData(sim)$Batch <- sim.cd[[match.idx[1]]]
        }
    }

    if ("Sample" %in% colnames(SummarizedExperiment::colData(sim))) {
        sample.lookup <- stats::setNames(prepared$sample.map$Sample, prepared$sample.map$vcf_id)
        sim.samples <- as.character(SummarizedExperiment::colData(sim)$Sample)
        mapped.samples <- unname(sample.lookup[sim.samples])
        keep.original <- is.na(mapped.samples)
        mapped.samples[keep.original] <- sim.samples[keep.original]
        SummarizedExperiment::colData(sim)$Sample <- mapped.samples
    }

    comparison.real.sce <- NULL
    comparison.sim <- NULL
    comparison.samples <- NULL
    if (!is.null(comparison.batch)) {
        keep.compare <- as.character(prepared$sce$Batch) == comparison.batch
        if (!is.null(comparison.pathology) &&
                "Pathology" %in% colnames(SummarizedExperiment::colData(prepared$sce))) {
            keep.compare <- keep.compare &
                as.character(prepared$sce$Pathology) == comparison.pathology
        }
        comparison.real.sce <- prepared$sce[, keep.compare, drop = FALSE]
        if (ncol(comparison.real.sce) == 0L) {
            stop(
                "No real-data cells remain for comparison.batch = '",
                comparison.batch, "'."
            )
        }

        comparison.samples <- unique(as.character(comparison.real.sce$Sample))
        comparison.aggregated <- simPICaggregateMeanByGroup(
            comparison.real.sce,
            comparison.real.sce$sample_batch
        )
        comparison.unit.sizes <- sort(table(comparison.real.sce$sample_batch), decreasing = TRUE)
        comparison.big.unit <- names(comparison.unit.sizes)[1]
        comparison.bigcounts <- comparison.real.sce[
            ,
            comparison.real.sce$sample_batch == comparison.big.unit,
            drop = FALSE
        ]
        comparison.sample.map <- prepared$sample.map[comparison.samples, , drop = FALSE]
        comparison.vcf <- prepared$vcf[, comparison.sample.map$vcf_id, drop = FALSE]
        comparison.ncells <- round(mean(as.integer(table(comparison.real.sce$Sample))))

        if (verbose) {
            message(
                "Simulating comparison subset for ", comparison.batch,
                if (!is.null(comparison.pathology)) paste0(" / ", comparison.pathology) else "",
                "..."
            )
        }

        comparison.params <- tryCatch(
            splatPopEstimatePeak(
                counts = comparison.bigcounts,
                means = getCounts(comparison.aggregated),
                params = params,
                sample.col = "Sample",
                batch.col = "Batch",
                aggregate.by = "sample_batch",
                pop.cv.bins = pop.cv.bins,
                verbose = FALSE
            ),
            error = function(e) {
                if (!grepl("negative 'bcv.common'", conditionMessage(e), fixed = TRUE)) {
                    stop(e)
                }
                splatPopEstimatePeak(
                    counts = comparison.bigcounts,
                    means = getCounts(comparison.aggregated),
                    params = params,
                    sample.col = "Sample",
                    batch.col = "Batch",
                    aggregate.by = "sample_batch",
                    pop.cv.bins = pop.cv.bins,
                    apply.bcv.correction = FALSE,
                    verbose = FALSE
                )
            }
        )
        comparison.params <- splatter::setParams(comparison.params, pop.quant.norm = FALSE)
        comparison.params <- splatter::setParams(
            comparison.params,
            batchCells = c(comparison.ncells),
            batch.size = 1L,
            eqtl.n = eqtl.n,
            similarity.scale = similarity.scale,
            eqtl.dist = eqtl.dist,
            ...
        )
        comparison.sim <- withCallingHandlers(
            splatPopSimulatePeak(
                params = comparison.params,
                vcf = comparison.vcf,
                gff = prepared$peak.annot,
                sparsify = sparsify,
                verbose = FALSE
            ),
            warning = function(w) {
                if (grepl("No SNP within eqtl.dist limit", conditionMessage(w), fixed = TRUE)) {
                    invokeRestart("muffleWarning")
                }
            }
        )
        comparison.cd <- SummarizedExperiment::colData(comparison.sim)
        if (!("Sample" %in% colnames(comparison.cd))) {
            sample.candidates <- c("Individual", "sample", "ind", "Donor")
            match.idx <- sample.candidates[sample.candidates %in% colnames(comparison.cd)]
            if (length(match.idx) > 0) {
                SummarizedExperiment::colData(comparison.sim)$Sample <- comparison.cd[[match.idx[1]]]
            }
        }
        if ("Sample" %in% colnames(SummarizedExperiment::colData(comparison.sim))) {
            sample.lookup <- stats::setNames(comparison.sample.map$Sample, comparison.sample.map$vcf_id)
            sim.samples <- as.character(SummarizedExperiment::colData(comparison.sim)$Sample)
            mapped.samples <- unname(sample.lookup[sim.samples])
            keep.original <- is.na(mapped.samples)
            mapped.samples[keep.original] <- sim.samples[keep.original]
            SummarizedExperiment::colData(comparison.sim)$Sample <- mapped.samples
        }
        SummarizedExperiment::colData(comparison.sim)$Batch <- comparison.batch
        if (!is.null(comparison.pathology)) {
            SummarizedExperiment::colData(comparison.sim)$Pathology <- comparison.pathology
        }
    }

    if (!is.null(comparison.real.sce) && !is.null(comparison.sim)) {
        real.plot.sce <- comparison.real.sce
        sim.plot.sce <- comparison.sim
        plot.samples <- comparison.samples
    } else {
        real.plot.sce <- simPICsubsetPcaSamples(
            prepared$sce,
            sample.col = "Sample",
            plot.samples = plot.samples,
            plot.n.samples = plot.n.samples
        )
        plot.samples <- unique(as.character(SummarizedExperiment::colData(real.plot.sce)$Sample))

        if ("Sample" %in% colnames(SummarizedExperiment::colData(sim))) {
            sim.plot.sce <- sim[, SummarizedExperiment::colData(sim)$Sample %in% plot.samples, drop = FALSE]
        } else {
            sim.plot.sce <- sim
        }
    }

    if (ncol(sim.plot.sce) == 0L) {
        stop("No simulated cells remain after subsetting to the plotting libraries.")
    }

    real.pca <- simPICplotPopulationPCA(
        real.plot.sce,
        sample.col = "Sample",
        batch.col = "Batch",
        pca.ntop = pca.ntop,
        pca.components = pca.components,
        verbose = FALSE
    )
    sim.pca <- simPICplotPopulationPCA(
        sim.plot.sce,
        sample.col = if ("Sample" %in% colnames(SummarizedExperiment::colData(sim.plot.sce))) "Sample" else NULL,
        batch.col = if ("Batch" %in% colnames(SummarizedExperiment::colData(sim.plot.sce))) "Batch" else NULL,
        pca.ntop = pca.ntop,
        pca.components = pca.components,
        verbose = FALSE
    )
    comparison.plots <- simPICplotPopulationComparison(
        real.pca = real.pca,
        sim.pca = sim.pca,
        title = sprintf(
            "Real and simulated cell-level PCA for the top %d libraries",
            length(plot.samples)
        )
    )
    bluster.comparison <- simPICplotBlusterComparison(
        real.sce = real.plot.sce,
        simulated.sce = sim.plot.sce,
        sample.col = "Sample",
        batch.col = if ("Batch" %in% colnames(SummarizedExperiment::colData(real.plot.sce))) "Batch" else NULL,
        subset.batch = comparison.batch,
        pathology.col = if ("Pathology" %in% colnames(SummarizedExperiment::colData(real.plot.sce))) "Pathology" else NULL,
        subset.pathology = comparison.pathology,
        plot.samples = if (is.null(comparison.batch)) plot.samples else NULL,
        plot.n.samples = plot.n.samples,
        pca.ntop = pca.ntop,
        pca.components = pca.components,
        point.size = 0.8,
        verbose = FALSE
    )

    return(list(
        sce = prepared$sce,
        bigcounts = prepared$bigcounts,
        aggregated = prepared$aggregated,
        params = params,
        sim = sim,
        comparison_real = comparison.real.sce,
        comparison_sim = comparison.sim,
        plot_samples = plot.samples,
        plots = list(
            real = real.pca$plots,
            simulated = sim.pca$plots,
            comparison = comparison.plots,
            bluster = bluster.comparison$plots
        ),
        sample_map = prepared$sample.map,
        peak_annot = prepared$peak.annot,
        vcf = prepared$vcf
    ))
}
