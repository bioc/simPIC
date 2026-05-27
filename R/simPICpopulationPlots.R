#' Plot Population-Style PCA Summaries for Real Data
#'
#' Create splatPop-style PCA plots from a real peak-by-cell matrix or
#' \code{SingleCellExperiment}. When metadata is not already available, sample
#' IDs can be inferred from cell names, which is useful for peak-by-cell
#' matrices where columns are encoded like \code{sample#barcode}.
#'
#' @param counts A peak-by-cell count matrix, sparse \code{Matrix}, or
#' \code{SingleCellExperiment}.
#' @param sample.col Optional name of a sample column already present in
#' \code{colData(counts)}. If \code{NULL}, sample IDs are inferred from column
#' names using \code{sample.pattern}.
#' @param batch.col Optional name of a batch column already present in
#' \code{colData(counts)}.
#' @param sample.pattern Regular expression used to strip the barcode suffix
#' from cell names when inferring sample IDs.
#' @param sample.replacement Replacement string used with
#' \code{sample.pattern}.
#' @param pca.ntop Number of most variable peaks to use for PCA.
#' @param pca.components Number of principal components to compute.
#' @param aggregate.by.sample Logical. Whether to also aggregate cells by sample
#' and create a donor-level PCA plot.
#' @param point.size Point size passed to \code{scater::plotPCA()}.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @return A list containing:
#' \describe{
#' \item{\code{cell_sce}}{Cell-level \code{SingleCellExperiment} with inferred
#' metadata, log-normalized counts, and PCA.}
#' \item{\code{sample_sce}}{Sample-aggregated \code{SingleCellExperiment} if
#' \code{aggregate.by.sample = TRUE}, otherwise \code{NULL}.}
#' \item{\code{plots}}{A named list of \code{ggplot} objects.}
#' }
#'
#' @examples
#' if (requireNamespace("scater", quietly = TRUE)) {
#'     counts <- matrix(rpois(50 * 24, lambda = 3), nrow = 50, ncol = 24)
#'     colnames(counts) <- paste0(
#'         rep(paste0("Sample", 1:6), each = 4),
#'         "#Cell",
#'         seq_len(ncol(counts))
#'     )
#'     plots <- simPICplotPopulationPCA(counts, verbose = FALSE)
#'     names(plots$plots)
#' }
#' @export
simPICplotPopulationPCA <- function(counts,
                                   sample.col = NULL,
                                   batch.col = NULL,
                                   sample.pattern = "#.*$",
                                   sample.replacement = "",
                                   pca.ntop = 2000,
                                   pca.components = 10,
                                   aggregate.by.sample = TRUE,
                                   point.size = 0.8,
                                   verbose = TRUE) {
    simPICrequireNamespace("scater", "simPICplotPopulationPCA()")
    simPICrequireNamespace("ggplot2", "simPICplotPopulationPCA()")

    checkmate::assertFlag(verbose)
    checkmate::assertFlag(aggregate.by.sample)
    checkmate::assertCount(pca.ntop, positive = TRUE)
    checkmate::assertCount(pca.components, positive = TRUE)
    checkmate::assertNumber(point.size, lower = 0, finite = TRUE)
    checkmate::assertString(sample.col, null.ok = TRUE, min.chars = 1)
    checkmate::assertString(batch.col, null.ok = TRUE, min.chars = 1)

    sce <- simPICcoerceCountsToSCE(counts)
    sce <- simPICfilterZeroLibraryCells(sce)

    col.names <- colnames(sce)
    if (is.null(col.names)) {
        col.names <- paste0("Cell", seq_len(ncol(sce)))
        colnames(sce) <- col.names
    }

    if (is.null(sample.col) || !(sample.col %in% colnames(SummarizedExperiment::colData(sce)))) {
        inferred <- simPICinferSampleFromColnames(
            col.names,
            pattern = sample.pattern,
            replacement = sample.replacement
        )
        SummarizedExperiment::colData(sce)$Sample <- inferred
        sample.col <- "Sample"
    }

    SummarizedExperiment::colData(sce)$Cell <- colnames(sce)
    SummarizedExperiment::colData(sce)$LibrarySize <- Matrix::colSums(getCounts(sce))

    if (verbose) {
        message("Computing log-normalized counts and PCA for cells...")
    }
    sce <- simPICaddLogcounts(sce)
    sce <- scater::runPCA(
        sce,
        exprs_values = "logcounts",
        ntop = min(pca.ntop, nrow(sce)),
        ncomponents = pca.components
    )
    cell.percent.var <- attr(SingleCellExperiment::reducedDim(sce, "PCA"), "percentVar")
    cell.title <- sprintf(
        "Cell-level PCA coloured by sample (PC1 %.2f%%, PC2 %.2f%%)",
        cell.percent.var[1],
        cell.percent.var[2]
    )

    plots <- list(
        cell_pca_by_sample = scater::plotPCA(
            sce,
            colour_by = sample.col,
            point_size = point.size
        ) +
            ggplot2::ggtitle(cell.title),
        cells_per_sample = ggplot2::ggplot(
            data.frame(Sample = SummarizedExperiment::colData(sce)[[sample.col]]),
            ggplot2::aes(x = .data$Sample)
        ) +
            ggplot2::geom_bar(fill = "#2C7FB8") +
            ggplot2::coord_flip() +
            ggplot2::labs(
                title = "Cells per sample",
                x = "Sample",
                y = "Cells"
            ) +
            ggplot2::theme_minimal(base_size = 11)
    )

    if (!is.null(batch.col) && batch.col %in% colnames(SummarizedExperiment::colData(sce))) {
        batch.levels <- unique(SummarizedExperiment::colData(sce)[[batch.col]])
        batch.levels <- batch.levels[!is.na(batch.levels)]

        if (length(batch.levels) <= 10L) {
            plots$cell_pca_by_batch <- scater::plotPCA(
                sce,
                colour_by = sample.col,
                shape_by = batch.col,
                point_size = point.size
            ) +
                ggplot2::ggtitle(
                    sprintf(
                        "Cell-level PCA coloured by sample and shaped by batch (PC1 %.2f%%, PC2 %.2f%%)",
                        cell.percent.var[1],
                        cell.percent.var[2]
                    )
                )
        }
    }

    sample.sce <- NULL
    if (aggregate.by.sample) {
        if (verbose) {
            message("Aggregating cells by sample and computing sample-level PCA...")
        }
        sample.sce <- simPICaggregateMeanByGroup(
            sce,
            SummarizedExperiment::colData(sce)[[sample.col]]
        )
        SummarizedExperiment::colData(sample.sce)$Sample <- colnames(sample.sce)
        sample.size.factors <- Matrix::colSums(getCounts(sample.sce))
        keep.samples <- sample.size.factors > 0

        if (!any(keep.samples)) {
            stop("No samples remain after aggregation with positive library size.")
        }

        sample.sce <- sample.sce[, keep.samples, drop = FALSE]
        sample.size.factors <- sample.size.factors[keep.samples]
        sample.sce <- simPICaddLogcounts(sample.sce,
            size.factors = sample.size.factors
        )
        sample.ncomp <- min(pca.components, max(2L, ncol(sample.sce) - 1L))
        sample.sce <- scater::runPCA(
            sample.sce,
            exprs_values = "logcounts",
            ntop = min(pca.ntop, nrow(sample.sce)),
            ncomponents = sample.ncomp
        )
        sample.percent.var <- attr(
            SingleCellExperiment::reducedDim(sample.sce, "PCA"),
            "percentVar"
        )

        plots$sample_pca <- scater::plotPCA(
            sample.sce,
            colour_by = "Sample",
            point_size = max(2, point.size * 2)
        ) +
            ggplot2::ggtitle(
                sprintf(
                    "Sample-level PCA from aggregated peak means (PC1 %.2f%%, PC2 %.2f%%)",
                    sample.percent.var[1],
                    sample.percent.var[2]
                )
            )
    }

    return(list(
        cell_sce = sce,
        sample_sce = sample.sce,
        plots = plots
    ))
}

simPICsubsetPcaSamples <- function(sce,
                                   sample.col = "Sample",
                                   plot.samples = NULL,
                                   plot.n.samples = 5L) {
    checkmate::assertClass(sce, "SingleCellExperiment")
    checkmate::assertString(sample.col, min.chars = 1)
    checkmate::assertCount(plot.n.samples, positive = TRUE)

    if (!(sample.col %in% colnames(SummarizedExperiment::colData(sce)))) {
        stop("'", sample.col, "' must exist in colData(sce).")
    }

    sample.ids <- as.character(SummarizedExperiment::colData(sce)[[sample.col]])
    if (is.null(plot.samples)) {
        sample.freq <- sort(table(sample.ids), decreasing = TRUE)
        plot.samples <- names(sample.freq)[seq_len(min(plot.n.samples, length(sample.freq)))]
    } else {
        plot.samples <- intersect(as.character(plot.samples), unique(sample.ids))
    }

    if (length(plot.samples) == 0L) {
        stop("No plot samples remain after subsetting.")
    }

    sce[, sample.ids %in% plot.samples, drop = FALSE]
}

simPICextractPcaFrame <- function(sce,
                                  sample.col = "Sample",
                                  dataset.label = "Real") {
    coords <- as.data.frame(SingleCellExperiment::reducedDim(sce, "PCA"))
    coords <- coords[, seq_len(min(2L, ncol(coords))), drop = FALSE]
    if (ncol(coords) < 2L) {
        stop("PCA output must contain at least two components.")
    }

    colnames(coords)[1:2] <- c("PC1", "PC2")
    coords$Sample <- as.character(SummarizedExperiment::colData(sce)[[sample.col]])
    coords$Dataset <- dataset.label
    coords
}

simPICplotPopulationComparison <- function(real.pca,
                                           sim.pca,
                                           title,
                                           point.size = 0.8) {
    simPICrequireNamespace("ggplot2", "simPICplotPopulationComparison()")

    checkmate::assertList(real.pca)
    checkmate::assertList(sim.pca)
    checkmate::assertString(title, min.chars = 1)
    checkmate::assertNumber(point.size, lower = 0, finite = TRUE)

    cell.df <- rbind(
        simPICextractPcaFrame(real.pca$cell_sce, dataset.label = "Real"),
        simPICextractPcaFrame(sim.pca$cell_sce, dataset.label = "Simulated")
    )

    plots <- list(
        cell_pca_by_sample = ggplot2::ggplot(
            cell.df,
            ggplot2::aes(x = .data$PC1, y = .data$PC2, colour = .data$Sample)
        ) +
            ggplot2::geom_point(alpha = 0.8, size = point.size) +
            ggplot2::facet_wrap(~Dataset, scales = "free") +
            ggplot2::labs(
                title = title,
                x = "PC1",
                y = "PC2",
                colour = "Sample"
            ) +
            ggplot2::theme_minimal(base_size = 11)
    )

    if (!is.null(real.pca$sample_sce) && !is.null(sim.pca$sample_sce)) {
        sample.df <- rbind(
            simPICextractPcaFrame(real.pca$sample_sce, dataset.label = "Real"),
            simPICextractPcaFrame(sim.pca$sample_sce, dataset.label = "Simulated")
        )

        plots$sample_pca <- ggplot2::ggplot(
            sample.df,
            ggplot2::aes(x = .data$PC1, y = .data$PC2, colour = .data$Sample)
        ) +
            ggplot2::geom_point(size = max(2, point.size * 2)) +
            ggplot2::facet_wrap(~Dataset, scales = "free") +
            ggplot2::labs(
                title = "Aggregated sample-level PCA for the same library subset",
                x = "PC1",
                y = "PC2",
                colour = "Sample"
            ) +
            ggplot2::theme_minimal(base_size = 11)
    }

    plots
}

#' Plot Poptrial-Style PCA, Silhouette, and Purity Comparisons
#'
#' Create side-by-side real-versus-simulated comparison panels inspired by the
#' bluster-based plotting workflow used in \code{Poptrial.Rmd}. The comparison
#' can be restricted either to a chosen set of samples or to a real-data batch,
#' in which case the simulated data are matched by the same sample IDs.
#'
#' @param real.sce Real-data \code{SingleCellExperiment}.
#' @param simulated.sce Simulated \code{SingleCellExperiment}.
#' @param sample.col Column in \code{colData()} containing sample IDs.
#' @param batch.col Optional column in the real-data \code{colData()} used to
#' define a batch or library subset such as \code{"Library5"}.
#' @param subset.batch Optional batch/library value used to subset the real data
#' before matching the simulated samples.
#' @param pathology.col Optional pathology column in the real-data
#' \code{colData()}.
#' @param subset.pathology Optional pathology value used together with
#' \code{subset.batch} to restrict the real-data comparison subset.
#' @param plot.samples Optional character vector of samples to compare. Ignored
#' when \code{subset.batch} is supplied.
#' @param plot.n.samples Number of top samples to keep when neither
#' \code{subset.batch} nor \code{plot.samples} is supplied.
#' @param pca.ntop Number of most variable peaks used for PCA.
#' @param pca.components Number of principal components to compute.
#' @param point.size Point size used in PCA panels.
#' @param verbose Logical. Whether to print progress messages.
#'
#' @return A list with the subsetted real and simulated SCEs, the selected
#' sample IDs, the bluster metric tables, and three ggplot panels:
#' \code{cell_pca_by_sample}, \code{silhouette_width}, and
#' \code{neighborhood_purity}.
#' @export
simPICplotBlusterComparison <- function(real.sce,
                                        simulated.sce,
                                        sample.col = "Sample",
                                        batch.col = NULL,
                                        subset.batch = NULL,
                                        pathology.col = NULL,
                                        subset.pathology = NULL,
                                        plot.samples = NULL,
                                        plot.n.samples = 5L,
                                        pca.ntop = 2000L,
                                        pca.components = 10L,
                                        point.size = 0.8,
                                        verbose = TRUE) {
    simPICrequireNamespace("bluster", "simPICplotBlusterComparison()")
    simPICrequireNamespace("ggplot2", "simPICplotBlusterComparison()")
    simPICrequireNamespace("scater", "simPICplotBlusterComparison()")

    checkmate::assertClass(real.sce, "SingleCellExperiment")
    checkmate::assertClass(simulated.sce, "SingleCellExperiment")
    checkmate::assertString(sample.col, min.chars = 1)
    checkmate::assertString(batch.col, null.ok = TRUE, min.chars = 1)
    checkmate::assertString(subset.batch, null.ok = TRUE, min.chars = 1)
    checkmate::assertString(pathology.col, null.ok = TRUE, min.chars = 1)
    checkmate::assertString(subset.pathology, null.ok = TRUE, min.chars = 1)
    checkmate::assertCharacter(plot.samples, null.ok = TRUE, min.chars = 1)
    checkmate::assertCount(plot.n.samples, positive = TRUE)
    checkmate::assertCount(pca.ntop, positive = TRUE)
    checkmate::assertCount(pca.components, positive = TRUE)
    checkmate::assertNumber(point.size, lower = 0, finite = TRUE)
    checkmate::assertFlag(verbose)

    if (!(sample.col %in% colnames(SummarizedExperiment::colData(real.sce)))) {
        stop("'", sample.col, "' must exist in the real-data colData.")
    }
    if (!(sample.col %in% colnames(SummarizedExperiment::colData(simulated.sce)))) {
        stop("'", sample.col, "' must exist in the simulated-data colData.")
    }

    if (!is.null(subset.batch)) {
        if (is.null(batch.col) ||
                !(batch.col %in% colnames(SummarizedExperiment::colData(real.sce)))) {
            stop(
                "A valid 'batch.col' in the real-data object is required ",
                "when 'subset.batch' is supplied."
            )
        }
        keep.real <- as.character(real.sce[[batch.col]]) == subset.batch
        if (!is.null(subset.pathology)) {
            if (is.null(pathology.col) ||
                    !(pathology.col %in% colnames(SummarizedExperiment::colData(real.sce)))) {
                stop(
                    "A valid 'pathology.col' in the real-data object is required ",
                    "when 'subset.pathology' is supplied."
                )
            }
            keep.real <- keep.real & as.character(real.sce[[pathology.col]]) == subset.pathology
        }
        real.sce <- real.sce[, keep.real, drop = FALSE]
        if (ncol(real.sce) == 0L) {
            stop("No real-data cells remain after subsetting to batch/library '", subset.batch, "'.")
        }
        plot.samples <- unique(as.character(real.sce[[sample.col]]))
    } else {
        real.sce <- simPICsubsetPcaSamples(
            real.sce,
            sample.col = sample.col,
            plot.samples = plot.samples,
            plot.n.samples = plot.n.samples
        )
        plot.samples <- unique(as.character(real.sce[[sample.col]]))
    }

    simulated.sce <- simulated.sce[, as.character(simulated.sce[[sample.col]]) %in% plot.samples, drop = FALSE]
    if (ncol(simulated.sce) == 0L) {
        stop("No simulated cells remain after matching the selected real-data samples.")
    }

    if (verbose) {
        message("Computing Poptrial-style PCA and bluster summaries...")
    }

    real.pca <- simPICaddLogcounts(simPICfilterZeroLibraryCells(real.sce))
    sim.pca <- simPICaddLogcounts(simPICfilterZeroLibraryCells(simulated.sce))
    real.pca <- scater::runPCA(
        real.pca,
        exprs_values = "logcounts",
        ntop = min(pca.ntop, nrow(real.pca)),
        ncomponents = pca.components
    )
    sim.pca <- scater::runPCA(
        sim.pca,
        exprs_values = "logcounts",
        ntop = min(pca.ntop, nrow(sim.pca)),
        ncomponents = pca.components
    )

    pca.df <- rbind(
        simPICextractPcaFrame(real.pca, sample.col = sample.col, dataset.label = "Real"),
        simPICextractPcaFrame(sim.pca, sample.col = sample.col, dataset.label = "Simulated")
    )

    kelly.colors <- c(
        "#FB9A99", "#1F78B4", "#FDBF6F", "#E31A1C", "#33A02C",
        "#FF7F00", "#6A3D9A", "#B15928", "#A6CEE3", "#B2DF8A",
        "#CAB2D6", "#191919", "#00C5CD", "#7FFF00", "#FF1493",
        "#FFD700", "#0000FF", "#8B4513", "#006400", "#4682B4"
    )
    real.samples <- sort(unique(as.character(real.pca[[sample.col]])))
    sim.samples <- sort(unique(as.character(sim.pca[[sample.col]])))
    n.colors <- max(length(real.samples), length(sim.samples))
    if (n.colors > length(kelly.colors)) {
        extra.colors <- grDevices::hcl.colors(n.colors - length(kelly.colors), "Dynamic")
        kelly.colors <- c(kelly.colors, extra.colors)
    }
    real.color.map <- stats::setNames(kelly.colors[seq_along(real.samples)], real.samples)
    sim.color.map <- stats::setNames(kelly.colors[seq_along(sim.samples)], sim.samples)
    library.color.map <- c(real.color.map, sim.color.map)

    sil.real <- bluster::approxSilhouette(
        SingleCellExperiment::reducedDim(real.pca, "PCA"),
        clusters = as.character(real.pca[[sample.col]])
    )
    sil.sim <- bluster::approxSilhouette(
        SingleCellExperiment::reducedDim(sim.pca, "PCA"),
        clusters = as.character(sim.pca[[sample.col]])
    )

    purity.real <- bluster::neighborPurity(
        SingleCellExperiment::reducedDim(real.pca, "PCA"),
        clusters = as.character(real.pca[[sample.col]])
    )
    purity.sim <- bluster::neighborPurity(
        SingleCellExperiment::reducedDim(sim.pca, "PCA"),
        clusters = as.character(sim.pca[[sample.col]])
    )

    sil.df <- rbind(
        data.frame(
            Width = sil.real$width,
            Cluster = as.factor(as.character(real.pca[[sample.col]])),
            Dataset = "Real"
        ),
        data.frame(
            Width = sil.sim$width,
            Cluster = as.factor(as.character(sim.pca[[sample.col]])),
            Dataset = "Simulated"
        )
    )
    purity.df <- rbind(
        data.frame(
            Purity = purity.real$purity,
            Cluster = as.factor(as.character(real.pca[[sample.col]])),
            Dataset = "Real"
        ),
        data.frame(
            Purity = purity.sim$purity,
            Cluster = as.factor(as.character(sim.pca[[sample.col]])),
            Dataset = "Simulated"
        )
    )

    plots <- list(
        cell_pca_by_sample = ggplot2::ggplot(
            pca.df,
            ggplot2::aes(x = .data$PC1, y = .data$PC2, colour = .data$Sample)
        ) +
            ggplot2::geom_point(alpha = 0.8, size = point.size) +
            ggplot2::facet_wrap(~Dataset, scales = "free") +
            ggplot2::scale_colour_manual(values = library.color.map) +
            ggplot2::labs(
                title = "Real and simulated PCA for the selected library/sample subset",
                x = "PC1",
                y = "PC2",
                colour = "Sample"
            ) +
            ggplot2::theme_minimal(base_size = 11) +
            ggplot2::theme(
                legend.position = "right",
                panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.4)
            ),
        silhouette_width = ggplot2::ggplot(
            sil.df,
            ggplot2::aes(x = .data$Cluster, y = .data$Width, fill = .data$Cluster)
        ) +
            ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
            ggplot2::geom_jitter(width = 0.2, size = 0.6, alpha = 0.5) +
            ggplot2::scale_fill_manual(values = library.color.map) +
            ggplot2::geom_hline(yintercept = 0, linetype = "dashed", colour = "black", linewidth = 0.8) +
            ggplot2::facet_wrap(~Dataset) +
            ggplot2::coord_cartesian(ylim = c(-0.30, 0.35)) +
            ggplot2::labs(
                title = "Approximate silhouette widths",
                x = "Sample",
                y = "Silhouette width"
            ) +
            ggplot2::theme_minimal(base_size = 11) +
            ggplot2::theme(
                legend.position = "none",
                axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.4)
            ),
        neighborhood_purity = ggplot2::ggplot(
            purity.df,
            ggplot2::aes(x = .data$Cluster, y = .data$Purity, fill = .data$Cluster)
        ) +
            ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7, width = 0.5) +
            ggplot2::geom_jitter(width = 0.2, size = 0.6, alpha = 0.5) +
            ggplot2::scale_fill_manual(values = library.color.map) +
            ggplot2::facet_wrap(~Dataset) +
            ggplot2::coord_cartesian(ylim = c(0, 1)) +
            ggplot2::labs(
                title = "Neighborhood purity",
                x = "Sample",
                y = "Purity"
            ) +
            ggplot2::theme_minimal(base_size = 11) +
            ggplot2::theme(
                legend.position = "none",
                axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
                panel.border = ggplot2::element_rect(colour = "black", fill = NA, linewidth = 0.4)
            )
    )

    list(
        real_sce = real.pca,
        simulated_sce = sim.pca,
        plot_samples = plot.samples,
        metrics = list(
            silhouette = sil.df,
            purity = purity.df
        ),
        plots = plots
    )
}
