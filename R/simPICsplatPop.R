#' Estimate Peak-Based splatPop Parameters
#'
#' Estimate splatPop population parameters from peak-by-cell counts. This is a
#' peak-based wrapper around \code{splatter::splatPopEstimate()} that can
#' optionally aggregate donor-level peak means from a
#' \code{SingleCellExperiment} object and then apply the simPIC BCV correction.
#'
#' @param counts Either a \code{SingleCellExperiment} object or a peak-by-cell
#' count matrix. When \code{means} is \code{NULL}, \code{counts} must be a
#' \code{SingleCellExperiment} with sample-level metadata so donor means can be
#' aggregated automatically.
#' @param means Optional dense matrix of aggregated peak means, with peaks in
#' rows and donors/samples in columns.
#' @param eqtl Optional empirical eQTL table to pass through to
#' \code{splatter::splatPopEstimate()}.
#' @param params Optional \code{SplatPopParams} object. If \code{NULL}, a new
#' object is created with \code{pop.cv.bins}.
#' @param sample.col Name of the sample/donor column in \code{colData(counts)}
#' when \code{counts} is a \code{SingleCellExperiment}.
#' @param batch.col Optional batch column in \code{colData(counts)} used when
#' \code{aggregate.by = "sample_batch"}.
#' @param aggregate.by Whether automatically derived means should be aggregated
#' by sample or by sample-batch combinations.
#' @param min.cells Minimum number of cells required per aggregation unit when
#' deriving \code{means} from a \code{SingleCellExperiment}.
#' @param pop.cv.bins Number of CV bins to use when \code{params} is created
#' internally.
#' @param apply.bcv.correction Logical. If \code{TRUE}, apply the simPIC BCV
#' correction rule directly to the estimated \code{SplatPopParams}.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @return A \code{SplatPopParams} object.
#'
#' @examples
#' if (requireNamespace("splatter", quietly = TRUE) &&
#'     requireNamespace("VariantAnnotation", quietly = TRUE)) {
#'     set.seed(101)
#'     gene_means <- rgamma(60, shape = 2, rate = 0.4)
#'     cell_scales <- runif(48, min = 0.7, max = 1.4)
#'     counts <- vapply(
#'         cell_scales,
#'         function(scale) {
#'             rnbinom(60, mu = gene_means * scale, size = 0.2)
#'         },
#'         numeric(60)
#'     )
#'     sce <- SingleCellExperiment::SingleCellExperiment(
#'         assays = list(counts = counts),
#'         colData = data.frame(
#'             Sample = rep(paste0("S", 1:12), each = 4),
#'             Batch = rep(c("B1", "B2"), each = 24)
#'         )
#'     )
#'     params <- splatPopEstimatePeak(
#'         counts = sce,
#'         sample.col = "Sample",
#'         batch.col = "Batch",
#'         aggregate.by = "sample",
#'         min.cells = 2,
#'         pop.cv.bins = 10,
#'         verbose = FALSE
#'     )
#'     params
#' }
#' @export
splatPopEstimatePeak <- function(counts = NULL,
                                means = NULL,
                                eqtl = NULL,
                                params = NULL,
                                sample.col = "Sample",
                                batch.col = NULL,
                                aggregate.by = c("sample", "sample_batch"),
                                min.cells = 1,
                                pop.cv.bins = 50,
                                apply.bcv.correction = TRUE,
                                verbose = TRUE) {
    simPICrequireNamespace("splatter", "splatPopEstimatePeak()")

    checkmate::assertFlag(verbose)
    checkmate::assertFlag(apply.bcv.correction)
    checkmate::assertString(sample.col, min.chars = 1)
    checkmate::assertString(batch.col, null.ok = TRUE, min.chars = 1)
    checkmate::assertCount(min.cells, positive = TRUE)

    aggregate.by <- match.arg(aggregate.by)

    if (is.null(params)) {
        params <- splatter::newSplatPopParams(pop.cv.bins = pop.cv.bins)
    } else {
        checkmate::assertClass(params, "SplatPopParams")
    }

    counts.mat <- NULL
    if (!is.null(counts)) {
        if (methods::is(counts, "SingleCellExperiment")) {
            counts <- simPICfilterZeroLibraryCells(counts)
            counts.mat <- simPICcoerceDenseMatrix(getCounts(counts), "counts")

            if (is.null(means)) {
                if (verbose) {
                    message("Aggregating peak means for splatPop estimation...")
                }
                means <- simPICaggregatePeakMeans(
                    counts,
                    sample.col = sample.col,
                    batch.col = batch.col,
                    aggregate.by = aggregate.by,
                    min.cells = min.cells
                )
            }
        } else {
            counts.mat <- simPICcoerceDenseMatrix(counts, "counts")
        }
    }

    if (is.null(counts.mat) && is.null(means)) {
        stop("At least one of 'counts' or 'means' must be supplied.")
    }

    if (is.null(means) && is.null(counts.mat)) {
        stop("'means' can only be omitted when 'counts' is provided.")
    }

    if (is.null(means) && !methods::is(counts, "SingleCellExperiment")) {
        stop(
            "Automatic mean aggregation requires 'counts' to be a ",
            "SingleCellExperiment."
        )
    }

    means <- simPICcoerceDenseMatrix(means, "means")

    if (!is.null(counts.mat) && !is.null(means) && nrow(counts.mat) != nrow(means)) {
        stop(
            "'counts' and 'means' must have the same number of peaks/rows."
        )
    }

    if (verbose) {
        message("Estimating splatPop peak parameters...")
    }
    params <- tryCatch(
        splatter::splatPopEstimate(
            counts = counts.mat,
            means = means,
            eqtl = eqtl,
            params = params
        ),
        error = function(e) {
            warning(
                "Count-based splatPop estimation failed (",
                conditionMessage(e),
                "). Falling back to mean-based population estimation.",
                call. = FALSE
            )

            fallback <- splatter::splatPopEstimate(
                counts = NULL,
                means = means,
                eqtl = eqtl,
                params = params
            )

            if (!is.null(counts.mat)) {
                fallback <- splatter::setParams(
                    fallback,
                    nGenes = nrow(counts.mat),
                    batchCells = c(ncol(counts.mat))
                )
            }

            fallback
        }
    )

    if (apply.bcv.correction && !is.null(counts.mat)) {
        params <- simPICapplySplatPopBCVCorrection(
            params,
            counts.mat,
            verbose = verbose
        )
    }

    return(params)
}


#' Simulate Peak-Based Population Data with splatPop
#'
#' Simulate peak-based single-cell data using \code{splatter::splatPopSimulate}
#' while coercing peak annotations into the GFF-like structure that splatPop
#' expects internally.
#'
#' @param params Optional \code{SplatPopParams} object. If \code{NULL}, a new
#' object is created.
#' @param vcf Optional \code{VariantAnnotation} VCF object. If \code{NULL}, the
#' \code{splatter} mock VCF is used.
#' @param method Simulation mode passed through to \code{splatPopSimulate()}.
#' @param gff Peak annotation supplied either as \code{NULL}, a GFF-like
#' \code{data.frame}, or a \code{GRanges}. Peak annotations are coerced to a
#' GFF-like table and stored as gene-like features for the underlying splatPop
#' machinery.
#' @param eqtl Optional empirical eQTL table.
#' @param means Optional empirical population means matrix.
#' @param key Optional splatPop key.
#' @param counts.only Logical. Whether to keep only counts in the simulated
#' object.
#' @param sparsify Logical. Whether to sparsify the resulting assays.
#' @param verbose Logical. Whether to print progress messages.
#' @param ... Additional parameters passed to \code{splatter::splatPopSimulate}.
#'
#' @return A \code{SingleCellExperiment} object.
#'
#' @examples
#' if (requireNamespace("splatter", quietly = TRUE) &&
#'     requireNamespace("VariantAnnotation", quietly = TRUE)) {
#'     params <- splatter::newSplatPopParams(nGenes = 20)
#'     params <- splatter::setParams(params, batchCells = c(20))
#'     sim <- splatPopSimulatePeak(
#'         params = params,
#'         vcf = splatter::mockVCF(),
#'         gff = splatter::mockGFF()[seq_len(20), ],
#'         sparsify = FALSE,
#'         verbose = FALSE
#'     )
#'     sim
#' }
#' @export
splatPopSimulatePeak <- function(params = NULL,
                                vcf = NULL,
                                method = c("single", "groups", "paths"),
                                gff = NULL,
                                eqtl = NULL,
                                means = NULL,
                                key = NULL,
                                counts.only = FALSE,
                                sparsify = TRUE,
                                verbose = TRUE,
                                ...) {
    simPICrequireNamespace("splatter", "splatPopSimulatePeak()")
    simPICrequireNamespace("VariantAnnotation", "splatPopSimulatePeak()")

    method <- match.arg(method)
    checkmate::assertFlag(counts.only)
    checkmate::assertFlag(sparsify)
    checkmate::assertFlag(verbose)

    means <- simPICcoerceDenseMatrix(means, "means")

    if (is.null(params)) {
        n.genes <- if (!is.null(gff)) nrow(as.data.frame(gff)) else 50
        params <- splatter::newSplatPopParams(nGenes = n.genes)
    } else {
        checkmate::assertClass(params, "SplatPopParams")
    }

    if (is.null(vcf)) {
        vcf <- splatter::mockVCF()
    }
    checkmate::assertClass(vcf, "CollapsedVCF")

    if (!is.null(gff)) {
        gff.n <- nrow(as.data.frame(gff))
        if (splatter::getParam(params, "nGenes") != gff.n) {
            params <- splatter::setParams(params, nGenes = gff.n)
        }
        gff <- simPICcoercePeakGFF(gff, n.features = gff.n)
    }

    sim <- splatter::splatPopSimulate(
        params = params,
        vcf = vcf,
        method = method,
        gff = gff,
        eqtl = eqtl,
        means = means,
        key = key,
        counts.only = counts.only,
        sparsify = sparsify,
        verbose = verbose,
        ...
    )

    S4Vectors::metadata(sim)$simPIC.population.model <- "splatPopPeak"
    return(sim)
}
