#' Compare SingleCellExperiment objects
#'
#' Combine data from several SingleCellExperiment objects and produce
#' some basic plots comparing them.
#'
#' @param sces named list of SingleCellExperiment objects to combine and
#'       compare.
#' @param point.size size of points in scatter plots.
#' @param point.alpha opacity of points in scatter plots.
#' @param fits whether to include fits in scatter plots.
#' @param colours vector of colours to use for each dataset.
#'
#' @details The returned list has three items:
#'
#' \describe{
#'    \item{\code{RowData}}{Combined row data from the provided
#'    SingleCellExperiments.}
#'    \item{\code{ColData}}{Combined column data from the provided
#'    SingleCellExperiments.}
#'    \item{\code{Plots}}{Comparison plots
#'        \describe{
#'            \item{\code{Means}}{Boxplot of mean distribution.}
#'            \item{\code{Variances}}{Boxplot of variance distribution.}
#'            \item{\code{MeanVar}}{Scatter plot with fitted lines showing
#'            the mean-variance relationship.}
#'            \item{\code{LibrarySizes}}{Boxplot of the library size
#'            distribution.}
#'            \item{\code{ZerosPeak}}{Boxplot of the percentage of each peak
#'            that is zero.}
#'            \item{\code{ZerosCell}}{Boxplot of the percentage of each cell
#'           that is zero.}
#'            \item{\code{MeanZeros}}{Scatter plot with fitted lines showing
#'            the mean-zeros relationship.}
#'    }
#'  }
#' }
#'
#' The plots returned by this function are created using
#' \code{\link[ggplot2]{ggplot}} and are only a sample of the kind of plots
#' you might like to consider. The data used to create these plots is also
#' returned and should be in the correct format to allow you to create
#' further plots using \code{\link[ggplot2]{ggplot}}.
#'
#' @return List containing the combined datasets and plots.
#' @examples
#' sim1 <- simPICsimulate(
#'     nPeaks = 1000, nCells = 500,
#'     pm.distr = "weibull", seed = 7856
#' )
#' sim2 <- simPICsimulate(
#'     nPeaks = 1000, nCells = 500,
#'     pm.distr = "gamma", seed = 4234
#' )
#' comparison <- simPICcompare(list(weibull = sim1, gamma = sim2))
#' names(comparison)
#' names(comparison$Plots)
#' @importFrom S4Vectors metadata<- metadata
#' @importFrom SingleCellExperiment cpm<- cpm
#' @importFrom SummarizedExperiment assay
#' @importFrom rlang .data
#' @importFrom stats quantile
#' @export
simPICcompare <- function(sces, point.size = 0.2, point.alpha = 0.1,
                        fits = TRUE, colours = NULL) {
    checkmate::assertList(sces,
        types = "SingleCellExperiment",
        any.missing = FALSE, min.len = 1, names = "unique"
    )
    checkmate::assertNumber(point.size, finite = TRUE)
    checkmate::assertNumber(point.alpha, lower = 0, upper = 1)
    checkmate::assertLogical(fits, any.missing = FALSE, len = 1)
    if (!is.null(colours)) {
        checkmate::assertCharacter(colours,
            any.missing = FALSE,
            len = length(sces)
        )
    } else {
        colours <- c(
            "#E41A1C", "#377EB8", "#4DAF4A", "#984EA3",
            "#008080", "#FFA07A", "#7B68EE", "#9ACD32", "black"
        )
    }

    for (name in names(sces)) {
        sce <- sces[[name]]
        rowData(sce)$Dataset <- name
        colData(sce)$Dataset <- name
        colData(sce)$sum <- colSums(assay(sce))
        sce <- scuttle::addPerCellQC(sce)
        sce <- scuttle::addPerFeatureQC(sce)
        sce <- addFeatureStats(sce, "counts")
        n.features <- colData(sce)$detected
        colData(sce)$PctZero <- 100 * (1 - n.features / nrow(sce))
        rowData(sce)$PctZero <- 100 - rowData(sce)$detected
        var.peaks <- rowData(sce)$var
        rowData(sce)$nzp <- 1 - (rowData(sce)$PctZero / 100)
        sces[[name]] <- sce
    }

    features <- rowData(sces[[1]])
    cells <- colData(sces[[1]])

    if (length(sces) > 1) {
        for (name in names(sces)[-1]) {
            sce <- sces[[name]]
            features <- rbindMatched(features, rowData(sce))
            cells <- rbindMatched(cells, colData(sce))
        }
    }

    features$Dataset <- factor(features$Dataset, levels = names(sces))
    cells$Dataset <- factor(cells$Dataset, levels = names(sces))
    features <- data.frame(features)
    cells <- data.frame(cells)

    means <- ggplot2::ggplot(
        features,
        ggplot2::aes(
            x = .data$Dataset, y = .data$mean,
            colour = .data$Dataset, fill = .data$Dataset
        )
    ) +
        ggplot2::geom_boxplot(
            width = 0.8, size = 0.6, alpha = 0.3,
            position = ggplot2::position_dodge2(0.5)
        ) +
        ggplot2::scale_colour_manual(values = colours) +
        ggplot2::scale_fill_manual(values = colours) +
        ggplot2::xlab("") +
        ggplot2::ylab("Peak means") +
        ggplot2::ggtitle("Distriution of peak means") +
        plot_theme()

    vars <- ggplot2::ggplot(
        features,
        ggplot2::aes(
            x = .data$Dataset, y = .data$VarCounts,
            colour = .data$Dataset, fill = .data$Dataset
        )
    ) +
        ggplot2::geom_boxplot(
            width = 0.8, size = 0.6, alpha = 0.3,
            position = ggplot2::position_dodge2(0.5)
        ) +
        ggplot2::scale_colour_manual(values = colours) +
        ggplot2::scale_fill_manual(values = colours) +
        ggplot2::xlab("") +
        ggplot2::ylab("Variance of peaks") +
        ggplot2::ggtitle("Distribution of peak variances") +
        plot_theme()

    mean.var <- ggplot2::ggplot(
        features,
        ggplot2::aes(
            x = .data$mean, y = .data$VarCounts,
            fill = .data$Dataset, color = .data$Dataset
        )
    ) +
        ggplot2::geom_point(size = point.size, alpha = point.alpha) +
        ggplot2::scale_colour_manual(values = colours) +
        ggplot2::scale_fill_manual(values = colours) +
        ggplot2::xlab("Peak-Means") +
        ggplot2::ylab("Peak Variance") +
        ggplot2::ggtitle("Mean-variance relationship") +
        plot_theme()

    libs <- ggplot2::ggplot(
        cells,
        ggplot2::aes(
            x = .data$Dataset, y = .data$sum, colour = .data$Dataset,
            fill = .data$Dataset
        )
    ) +
        ggplot2::geom_boxplot(
            width = 0.8, size = 0.6, alpha = 0.3,
            position = ggplot2::position_dodge2(0.5)
        ) +
        ggplot2::scale_y_continuous(labels = scales::comma) +
        ggplot2::scale_colour_manual(values = colours) +
        ggplot2::scale_fill_manual(values = colours) +
        ggplot2::xlab("") +
        ggplot2::ylab("library size") +
        ggplot2::ggtitle("Distribution of library sizes") +
        plot_theme()

    z.peak <- ggplot2::ggplot(
        features,
        ggplot2::aes(
            x = .data$Dataset, y = .data$PctZero,
            color = .data$Dataset, fill = .data$Dataset
        )
    ) +
        ggplot2::geom_boxplot(
            width = 0.8, size = 0.6, alpha = 0.3,
            position = ggplot2::position_dodge2(0.5)
        ) +
        ggplot2::scale_y_continuous(limits = c(0, 100)) +
        ggplot2::scale_colour_manual(values = colours) +
        ggplot2::scale_fill_manual(values = colours) +
        ggplot2::xlab("") +
        ggplot2::ylab("% zeros per peak") +
        ggplot2::ggtitle("Distribution of zeros per peak") +
        plot_theme()

    z.cell <- ggplot2::ggplot(
        cells,
        ggplot2::aes(
            x = .data$Dataset, y = .data$PctZero,
            colour = .data$Dataset, fill = .data$Dataset
        )
    ) +
        ggplot2::geom_boxplot(
            width = 0.8, size = 0.6, alpha = 0.3,
            position = ggplot2::position_dodge2(0.5)
        ) +
        ggplot2::ylim(50, quantile(features$PctZero, 1)) +
        ggplot2::scale_colour_manual(values = colours) +
        ggplot2::scale_fill_manual(values = colours) +
        ggplot2::xlab("") +
        ggplot2::ylab("% zeros per cell") +
        ggplot2::ggtitle("Distribution of zeros per cell") +
        plot_theme()

    mean.zeros <- ggplot2::ggplot(
        features,
        ggplot2::aes(
            x = .data$mean, y = .data$PctZero,
            color = .data$Dataset, fill = .data$Dataset
        )
    ) +
        ggplot2::geom_point(size = point.size, alpha = point.alpha) +
        ggplot2::scale_x_log10(labels = scales::comma) +
        ggplot2::scale_colour_manual(values = colours) +
        ggplot2::scale_fill_manual(values = colours) +
        ggplot2::xlab("Peak-Means") +
        ggplot2::ylab("% zeros") +
        ggplot2::ggtitle("Mean-zeros relationship") +
        plot_theme()

    pm.nzp <- ggplot2::ggplot(
        features,
        ggplot2::aes(
            x = .data$nzp, y = .data$mean,
            fill = .data$Dataset, color = .data$Dataset
        )
    ) +
        ggplot2::geom_point(size = point.size, alpha = point.alpha) +
        ggplot2::scale_colour_manual(values = colours) +
        ggplot2::scale_fill_manual(values = colours) +
        ggplot2::xlab("Non-zero proportion") +
        ggplot2::ylab("Peak-Means") +
        ggplot2::ggtitle("Correlation- peak mean and non-zero proportion") +
        plot_theme()

    if (fits) {
        mean.var <- mean.var +
            ggplot2::geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"))
        pm.nzp <- pm.nzp +
            ggplot2::geom_smooth(method = "gam", formula = y ~ s(x, bs = "cs"))
    }
    comparison <- list(
        RowData = features, ColData = cells,
        Plots = list(
            Means = means, Variances = vars, MeanVar = mean.var,
            LibrarySizes = libs, ZerosPeak = z.peak, ZerosCell = z.cell,
            MeanZeros = mean.zeros, PeakMeanNzp = pm.nzp
        )
    )
    return(comparison)
}

#' Add feature statistics
#'
#' Add additional feature statistics to a SingleCellExperiment object
#'
#' @param sce SingleCellExperiment to add feature statistics to.
#' @param value the count value to calculate statistics.
#' @param log logical. Whether to take log2 before calculating statistics.
#' @param offset offset to add to avoid taking log of zero.
#' @param no.zeros logical. Whether to remove all zeros from each feature
#'        before calculating statistics.
#'
#' @details
#' Currently adds the following statistics: mean and variance. Statistics
#' are added to the \code{\link{rowData}} slot and are named
#' \code{Stat[Log]Value[No0]} where \code{Log} and \code{No0} are added if
#' those arguments are true.
#' @return SingleCellExperiment with additional feature statistics
#' @importFrom SummarizedExperiment rowData rowData<-
#' @importFrom methods is
addFeatureStats <- function(sce, value = "counts",
                            log = FALSE, offset = 1, no.zeros = FALSE) {
    checkmate::assertClass(sce, "SingleCellExperiment")
    checkmate::assertLogical(log)
    checkmate::assertNumber(offset, lower = 0)
    checkmate::assertLogical(no.zeros)
    value <- match.arg(value)

    switch(value,
        counts = {
            values <- BiocGenerics::counts(sce)
            suffix <- "Counts"
        }
    )

    if (is(values, "dgCMatrix")) {
        values <- as.matrix(values)
    }

    if (no.zeros) {
        values[values == 0] <- NA
        suffix <- paste0(suffix, "No0")
    }

    if (log) {
        values <- log2(values + offset)
        suffix <- paste0("Log", suffix)
    }

    mean.str <- paste0("Mean", suffix)
    var.str <- paste0("Var", suffix)

    rowData(sce)[, mean.str] <- rowMeans(values, na.rm = TRUE)
    rowData(sce)[, var.str] <- matrixStats::rowVars(values, na.rm = TRUE)

    return(sce)
}

#' Custom theme for ggplot2
#'
#' This function defines a custom theme for ggplot2 to ensure consistent
#' visual appearance across multiple plots.
#' @return A ggplot2 theme object with predefined settings.
plot_theme <- function() {
    ggplot2::theme(
        axis.text = ggplot2::element_text(size = 12, colour = "black"),
        axis.title = ggplot2::element_text(size = 12),
        panel.grid = ggplot2::element_blank(),
        panel.border = ggplot2::element_rect(
            color = "black", fill = NA,
            size = 1
        ),
        panel.background = ggplot2::element_blank(),
        plot.background = ggplot2::element_blank(),
        plot.title = ggplot2::element_text(face = "bold", color = "black"),
        legend.position = "none"
    )
}
