#' Create a newsimPICcount object to store parameters.
#' 
#' This function creates the object variable which is passed in all
#' functions.
#'
#' @param ... Variables to set newsimPICcount object parameters.
#' @return new object from class simPICcount.
#'
#' @examples
#' object <- newsimPICcount()
#'
#' @importFrom methods new
#' @title newsimPICcount
#' @name newsimPICcount
#' @export
newsimPICcount <- function(...) {
    object <- methods::new("simPICcount")
    object <- setsimPICparameters(object, ...)

    return(object)
}

#' Set simPIC parameters
#'
#' Set input parameters of the simPICcount object.
#'
#' @param object input simPICcount object.
#' @param update new parameters.
#' @param ... set new parameters for simPICcount object.
#'
#' @return simPICcount object with updated parameters.
#'
#' @examples
#' object <- newsimPICcount()
#' object <- setsimPICparameters(object, nCells = 200, nPeaks = 500)
#'
#' @importFrom methods slot<- validObject
#' @export
setsimPICparameters <- function(object, update = NULL, ...) {
    checkmate::assertClass(object, classes = "simPICcount")
    checkmate::assertList(update, null.ok = TRUE)

    update <- c(update, list(...))

    if (length(update) > 0) {
        for (name in names(update)) {
            value <- update[[name]]
            checkmate::assertString(name)

            slot(object, name) <- value
            validObject(object)
        }
    }
    return(object)
}

#' Get a single simPICcount parameter
#'
#' Get the value of a single variable from input simPICcount object.
#'
#' @param object input simPICcount object.
#' @param name name of the parameter.
#'
#' @return Value of the input parameter.
#'
#' @examples
#' object <- newsimPICcount()
#' nPeaks <- simPICget(object, "nPeaks")
#'
#' @importFrom methods slot
#' @export
simPICget <- function(object, name) {
    slot(object, name)
}

#' Get parameters
#'
#' Get multiple parameter values from a simPIC object.
#'
#' @param object input object to get values from.
#' @param names vector of names of the parameters to get.
#'
#' @return List with the values of the selected parameters.
#' @examples
#' object <- newsimPICcount()
#' simPICgetparameters(object, c("nPeaks", "nCells", "peak.mean.shape"))
#' @export
simPICgetparameters <- function(object, names) {
    checkmate::assertClass(object, classes = "simPICcount")
    checkmate::assertCharacter(names, min.len = 1, any.missing = FALSE)

    params.list <- lapply(names, simPICget, object = object)
    names(params.list) <- names

    return(params.list)
}
#' Get counts from Single Cell Experiment object
#'
#' Get counts matrix from a SingleCellExperiment object. If counts is
#' missing a warning is issued and the first assay is returned.
#'
#' @param sce SingleCellExperiment object
#' @return counts matrix
getCounts <- function(sce) {
    checkmate::assertClass(sce, "SingleCellExperiment")

    if ("counts" %in% SummarizedExperiment::assayNames(sce)) {
        counts <- SingleCellExperiment::counts(sce)
    } else {
        warning("counts assay is missing, using the first assay instead")
        counts <- SummarizedExperiment::assay(sce)
    }

    return(counts)
}

#' Ensure counts assay is first
#'
#' Reorders assays in a SingleCellExperiment so that the \code{counts} assay
#' appears first. This keeps \code{assay(sce)} behavior predictable for
#' downstream code that still assumes the primary assay is counts.
#'
#' @param sce SingleCellExperiment object.
#' @return SingleCellExperiment with counts assay first if present.
ensureCountsFirst <- function(sce) {
    checkmate::assertClass(sce, "SingleCellExperiment")

    assay_names <- SummarizedExperiment::assayNames(sce)
    if (!"counts" %in% assay_names || identical(assay_names[[1]], "counts")) {
        return(sce)
    }

    new_order <- c("counts", setdiff(assay_names, "counts"))
    SummarizedExperiment::assays(sce, withDimnames = FALSE) <-
        SummarizedExperiment::assays(sce, withDimnames = FALSE)[new_order]

    sce
}

simPICaggregateMeanByGroup <- function(sce, ids) {
    checkmate::assertClass(sce, "SingleCellExperiment")
    checkmate::assertAtomicVector(ids, len = ncol(sce))

    keep <- !is.na(ids) & ids != ""
    if (!any(keep)) {
        stop("No cells remain after removing missing aggregation identifiers.")
    }

    sce <- sce[, keep, drop = FALSE]
    ids <- factor(ids[keep], levels = unique(ids[keep]))
    counts <- getCounts(sce)
    design <- Matrix::sparse.model.matrix(~ 0 + ids)
    colnames(design) <- levels(ids)

    summed <- counts %*% design
    group.sizes <- as.numeric(table(ids)[colnames(design)])
    means <- summed %*% Matrix::Diagonal(x = 1 / group.sizes)
    rownames(means) <- rownames(sce)
    colnames(means) <- colnames(design)

    SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = means),
        rowData = SummarizedExperiment::rowData(sce)
    )
}

simPICaddLogcounts <- function(sce, size.factors = NULL) {
    checkmate::assertClass(sce, "SingleCellExperiment")

    counts <- getCounts(sce)
    if (is.null(size.factors)) {
        lib.size <- Matrix::colSums(counts)
        positive <- lib.size > 0
        if (!all(positive)) {
            stop("All cells must have positive library sizes for normalization.")
        }
        size.factors <- lib.size / mean(lib.size)
    }

    checkmate::assertNumeric(size.factors,
        len = ncol(sce), lower = 0,
        any.missing = FALSE, finite = TRUE
    )
    if (any(size.factors == 0)) {
        stop("size.factors must be positive.")
    }

    logcounts <- log2(counts %*% Matrix::Diagonal(x = 1 / size.factors) + 1)
    dimnames(logcounts) <- dimnames(sce)
    SummarizedExperiment::assay(sce, "logcounts") <- logcounts
    sce
}

#' Bind rows (matched)
#'
#' Bind the rows of two data frames, keeping only the columns that are
#' common to both.
#'
#' @param df1 first data.frame to bind.
#' @param df2 second data.frame to bind.
#'
#' @return data.frame containing rows from \code{df1} and \code{df2} but
#' only common columns.
rbindMatched <- function(df1, df2) {
    common.names <- intersect(colnames(df1), colnames(df2))
    if (length(common.names) < 2) {
        stop("There must be at least two columns in common")
    }
    combined <- rbind(df1[, common.names], df2[, common.names])

    return(combined)
}

#' Convert Sparse Matrix to SingleCellExperiment object
#'
#' This function converts a sparse matrix into a SingleCellExperiment (SCE)
#' object.
#'
#' @param sparse_data A sparse matrix containing count data, where rows are
#' peaks and columns represent cells.
#' @return A SingleCellExperiment (SCE) object with the sparse matrix stored
#' in the "counts" assay.
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom Matrix Matrix
convert_to_SCE <- function(sparse_data) {
    sce <- SingleCellExperiment(assays = list(counts = sparse_data))
    rownames(sce) <- paste0("Peak", seq_len(nrow(sparse_data)))
    colnames(sce) <- paste0("Cell", seq_len(ncol(sparse_data)))
    return(sce)
}

simPICrequireNamespace <- function(pkg, fn) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
        stop(
            "Package '", pkg, "' is required for ", fn,
            ". Please install it first.",
            call. = FALSE
        )
    }
}

simPICcoerceDenseMatrix <- function(x, name) {
    if (is.null(x)) {
        return(NULL)
    }

    if (methods::is(x, "SingleCellExperiment")) {
        x <- getCounts(x)
    }

    if (is.data.frame(x)) {
        x <- as.matrix(x)
    }

    if (methods::is(x, "Matrix")) {
        warning(
            name, " is sparse; converting to a dense matrix for splatPop.",
            call. = FALSE
        )
        x <- as.matrix(x)
    }

    if (!is.matrix(x)) {
        stop(name, " must be a matrix, data.frame, or SingleCellExperiment.")
    }

    storage.mode(x) <- "double"

    if (anyNA(x) || any(!is.finite(x))) {
        stop(name, " must not contain missing or non-finite values.")
    }

    if (any(x < 0)) {
        stop(name, " must contain non-negative values.")
    }

    return(x)
}

simPICfilterZeroLibraryCells <- function(sce) {
    counts <- getCounts(sce)
    keep <- Matrix::colSums(counts) > 0

    if (!any(keep)) {
        stop("No cells remain after removing zero-library cells.")
    }

    sce <- sce[, keep, drop = FALSE]
    return(sce)
}

simPICaggregatePeakMeans <- function(sce,
                                    sample.col,
                                    batch.col = NULL,
                                    aggregate.by = c("sample", "sample_batch"),
                                    min.cells = 1) {
    simPICrequireNamespace("scuttle", "simPICaggregatePeakMeans()")

    aggregate.by <- match.arg(aggregate.by)

    cd <- SummarizedExperiment::colData(sce)
    cd.names <- colnames(cd)

    if (!(sample.col %in% cd.names)) {
        stop("'", sample.col, "' is not a column in colData(counts).")
    }

    ids <- as.character(cd[[sample.col]])

    if (aggregate.by == "sample_batch") {
        if (is.null(batch.col)) {
            stop("'batch.col' must be supplied when aggregate.by = 'sample_batch'.")
        }
        if (!(batch.col %in% cd.names)) {
            stop("'", batch.col, "' is not a column in colData(counts).")
        }
        ids <- paste(ids, as.character(cd[[batch.col]]), sep = "_")
    }

    ids[is.na(ids) | ids == ""] <- NA_character_
    keep.cells <- !is.na(ids)
    if (!any(keep.cells)) {
        stop("No cells remain after removing missing aggregation identifiers.")
    }

    sce <- sce[, keep.cells, drop = FALSE]
    ids <- ids[keep.cells]

    unit.sizes <- table(ids)
    keep.units <- names(unit.sizes)[unit.sizes >= min.cells]

    if (length(keep.units) == 0) {
        stop("No aggregation units remain after applying min.cells = ", min.cells, ".")
    }

    keep.cells <- ids %in% keep.units
    sce <- sce[, keep.cells, drop = FALSE]
    ids <- ids[keep.cells]

    aggregated <- simPICaggregateMeanByGroup(sce, ids)

    means <- simPICcoerceDenseMatrix(getCounts(aggregated), "means")
    return(means)
}

simPICapplySplatPopBCVCorrection <- function(params, counts, verbose = TRUE) {
    checkmate::assertClass(params, "SplatPopParams")
    counts <- simPICcoerceDenseMatrix(counts, "counts")

    design <- matrix(1, ncol(counts), 1)
    disps <- edgeR::estimateDisp(counts, design = design)
    bcv.common <- -0.3 + 0.15 * disps$common.dispersion

    if (bcv.common < 0) {
        stop(
            "The simPIC BCV correction produced a negative 'bcv.common' (",
            signif(bcv.common, 4),
            "), which is invalid for SplatPopParams. ",
            "Please use counts with higher dispersion or override the ",
            "BCV settings manually."
        )
    }

    params <- splatter::setParams(
        params,
        bcv.common = bcv.common,
        bcv.df = disps$prior.df
    )

    if (verbose) {
        message("Applied simPIC BCV correction to splatPop parameters.")
    }

    return(params)
}

simPICcoercePeakGFF <- function(gff, n.features = NULL) {
    if (is.null(gff)) {
        return(NULL)
    }

    if (methods::is(gff, "GRanges")) {
        gff <- as.data.frame(gff)
    } else {
        checkmate::assertDataFrame(gff, min.rows = 1)
        gff <- as.data.frame(gff, stringsAsFactors = FALSE)
    }

    seqnames <- starts <- ends <- ids <- NULL

    if (all(c("seqnames", "start", "end") %in% colnames(gff))) {
        seqnames <- gff[["seqnames"]]
        starts <- gff[["start"]]
        ends <- gff[["end"]]
        ids <- if ("peak_id" %in% colnames(gff)) gff[["peak_id"]] else rownames(gff)
    } else if (all(c("chromosome", "geneStart", "geneEnd") %in% colnames(gff))) {
        seqnames <- gff[["chromosome"]]
        starts <- gff[["geneStart"]]
        ends <- gff[["geneEnd"]]
        ids <- rownames(gff)
    } else if (all(c("V1", "V4", "V5") %in% colnames(gff))) {
        seqnames <- gff[["V1"]]
        starts <- gff[["V4"]]
        ends <- gff[["V5"]]
        ids <- if ("V9" %in% colnames(gff)) gff[["V9"]] else rownames(gff)
    } else if (ncol(gff) >= 5) {
        seqnames <- gff[[1]]
        starts <- gff[[4]]
        ends <- gff[[5]]
        ids <- if (ncol(gff) >= 9) gff[[9]] else rownames(gff)
    } else {
        stop(
            "Peak annotation must be a GRanges or data.frame with either ",
            "seqnames/start/end columns or GFF-like columns."
        )
    }

    if (is.null(ids) || all(is.na(ids)) || length(ids) == 0) {
        ids <- paste0("peak_", seq_len(nrow(gff)))
    }

    starts <- as.integer(as.character(starts))
    ends <- as.integer(as.character(ends))

    start.fixed <- pmin(starts, ends, na.rm = TRUE)
    end.fixed <- pmax(starts, ends, na.rm = TRUE)

    out <- data.frame(
        V1 = as.character(seqnames),
        V2 = "simPIC",
        V3 = "gene",
        V4 = start.fixed,
        V5 = end.fixed,
        V6 = ".",
        V7 = ".",
        V8 = ".",
        V9 = paste0("ID=", make.unique(as.character(ids))),
        stringsAsFactors = FALSE
    )

    if (anyNA(out$V4) || anyNA(out$V5)) {
        stop("Peak annotation start/end positions must be finite integers.")
    }

    if (!is.null(n.features) && nrow(out) != n.features) {
        stop(
            "Peak annotation has ", nrow(out), " rows but expected ",
            n.features, "."
        )
    }

    return(out)
}

simPICcoerceCountsToSCE <- function(counts, assay.name = "counts") {
    checkmate::assertString(assay.name, min.chars = 1)

    if (methods::is(counts, "SingleCellExperiment")) {
        return(counts)
    }

    if (is.data.frame(counts)) {
        counts <- as.matrix(counts)
    }

    if (!is.matrix(counts) && !methods::is(counts, "Matrix")) {
        stop(
            "'counts' must be a matrix, sparse Matrix, data.frame, or ",
            "SingleCellExperiment."
        )
    }

    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = stats::setNames(list(counts), assay.name)
    )

    if (!is.null(rownames(counts))) {
        rownames(sce) <- rownames(counts)
    }
    if (!is.null(colnames(counts))) {
        colnames(sce) <- colnames(counts)
    }

    return(sce)
}

simPICinferSampleFromColnames <- function(cells,
                                         pattern = "#.*$",
                                         replacement = "") {
    checkmate::assertCharacter(cells, any.missing = FALSE)
    checkmate::assertString(pattern, min.chars = 1)
    checkmate::assertString(replacement)

    samples <- sub(pattern, replacement, cells)

    if (anyNA(samples) || any(samples == "")) {
        stop(
            "Unable to infer sample IDs from cell names using pattern '",
            pattern, "'."
        )
    }

    return(samples)
}

simPICmicrogliaExamplePaths <- function() {
    extdir <- system.file("extdata", package = "simPIC")

    if (identical(extdir, "")) {
        extdir <- file.path(getwd(), "inst", "extdata")
    }

    list(
        sce = file.path(extdir, "microglia_sce.rds"),
        sample_map = file.path(extdir, "microglia_sample_map.tsv"),
        peak_annot = file.path(extdir, "microglia_peak_annot_chr14.tsv"),
        vcf = file.path(extdir, "microglia_chr14.vcf")
    )
}

simPICreadMicrogliaSampleMap <- function(sample.map) {
    if (is.character(sample.map) && length(sample.map) == 1L) {
        sample.map <- utils::read.delim(
            sample.map,
            header = TRUE,
            stringsAsFactors = FALSE
        )
    }

    checkmate::assertDataFrame(sample.map, min.rows = 1)
    sample.map <- as.data.frame(sample.map, stringsAsFactors = FALSE)

    if (ncol(sample.map) < 2) {
        stop("'sample.map' must contain at least two columns.")
    }

    colnames(sample.map)[seq_len(2)] <- c("vcf_id", "Sample")
    sample.map <- sample.map[, c("vcf_id", "Sample"), drop = FALSE]
    sample.map <- sample.map[!duplicated(sample.map$vcf_id), , drop = FALSE]
    sample.map <- sample.map[!is.na(sample.map$vcf_id) &
        !is.na(sample.map$Sample) &
        sample.map$vcf_id != "" &
        sample.map$Sample != "", , drop = FALSE]

    if (nrow(sample.map) == 0) {
        stop("No valid rows remain in 'sample.map'.")
    }

    rownames(sample.map) <- sample.map$Sample
    return(sample.map)
}

simPICreadPeakAnnotation <- function(peak.annot) {
    if (is.character(peak.annot) && length(peak.annot) == 1L) {
        peak.annot <- utils::read.delim(
            peak.annot,
            header = TRUE,
            stringsAsFactors = FALSE
        )
    }

    checkmate::assertDataFrame(peak.annot, min.rows = 1)
    peak.annot <- as.data.frame(peak.annot, stringsAsFactors = FALSE)

    if (!("peak_id" %in% colnames(peak.annot))) {
        peak.annot$peak_id <- paste0("Peak", seq_len(nrow(peak.annot)))
    }

    rownames(peak.annot) <- peak.annot$peak_id
    return(peak.annot)
}

simPICreadMicrogliaVcf <- function(vcf) {
    if (is.character(vcf) && length(vcf) == 1L) {
        simPICrequireNamespace("VariantAnnotation", "simPICreadMicrogliaVcf()")
        vcf <- VariantAnnotation::readVcf(vcf)
    }

    checkmate::assertClass(vcf, "CollapsedVCF")
    return(vcf)
}

simPICsanitizePeakAnnotationForRowData <- function(peak.annot) {
    peak.annot <- as.data.frame(peak.annot, stringsAsFactors = FALSE)
    reserved <- c(
        "seqnames", "ranges", "strand", "start",
        "end", "width", "element"
    )
    rename.idx <- colnames(peak.annot) %in% reserved
    colnames(peak.annot)[rename.idx] <- paste0(
        "peak_",
        colnames(peak.annot)[rename.idx]
    )
    peak.annot
}

simPICsubsetRowsByPeakAnnotation <- function(sce, peak.annot) {
    checkmate::assertClass(sce, "SingleCellExperiment")
    peak.annot <- simPICreadPeakAnnotation(peak.annot)

    keep <- rownames(sce) %in% peak.annot$peak_id

    if (!any(keep)) {
        stop("No overlapping peaks remain after applying the peak annotation.")
    }

    sce <- sce[keep, , drop = FALSE]
    peak.annot <- peak.annot[rownames(sce), , drop = FALSE]
    safe.annot <- simPICsanitizePeakAnnotationForRowData(peak.annot)
    SummarizedExperiment::rowData(sce) <- S4Vectors::DataFrame(safe.annot)

    return(list(sce = sce, peak.annot = peak.annot))
}

simPICalignVcfToSamples <- function(vcf, sample.map, samples) {
    sample.map <- simPICreadMicrogliaSampleMap(sample.map)
    keep.map <- sample.map[sample.map$Sample %in% samples, , drop = FALSE]

    if (nrow(keep.map) == 0) {
        stop("No mapped samples remain after filtering the Microglia example.")
    }

    vcf.samples <- colnames(VariantAnnotation::geno(vcf)$GT)
    keep.map <- keep.map[keep.map$vcf_id %in% vcf.samples, , drop = FALSE]

    if (nrow(keep.map) == 0) {
        stop("No VCF samples remain after applying the Microglia sample map.")
    }

    vcf <- vcf[, keep.map$vcf_id]
    rownames(keep.map) <- keep.map$Sample

    return(list(vcf = vcf, sample.map = keep.map))
}

simPICharmonizePeakSeqnames <- function(peak.annot, vcf) {
    peak.annot <- simPICreadPeakAnnotation(peak.annot)
    vcf.seqlevels <- as.character(
        GenomeInfoDb::seqlevels(SummarizedExperiment::rowRanges(vcf))
    )
    vcf.seqlevels <- vcf.seqlevels[!is.na(vcf.seqlevels) & vcf.seqlevels != ""]

    if (length(vcf.seqlevels) == 0 || !("seqnames" %in% colnames(peak.annot))) {
        return(peak.annot)
    }

    if (all(grepl("^chr", peak.annot$seqnames)) && !any(grepl("^chr", vcf.seqlevels))) {
        peak.annot$seqnames <- sub("^chr", "", peak.annot$seqnames)
    } else if (!any(grepl("^chr", peak.annot$seqnames)) && any(grepl("^chr", vcf.seqlevels))) {
        peak.annot$seqnames <- paste0("chr", peak.annot$seqnames)
    }

    rownames(peak.annot) <- peak.annot$peak_id
    peak.annot
}

simPICprepareMicrogliaExample <- function(sce,
                                         sample.map,
                                         peak.annot,
                                         vcf,
                                         min.cells = 100L,
                                         force.keep.batch = NULL,
                                         force.keep.pathology = NULL) {
    checkmate::assertClass(sce, "SingleCellExperiment")
    checkmate::assertCount(min.cells, positive = TRUE)
    checkmate::assertString(force.keep.batch, null.ok = TRUE, min.chars = 1)
    checkmate::assertString(force.keep.pathology, null.ok = TRUE, min.chars = 1)

    sce <- SingleCellExperiment::SingleCellExperiment(
        assays = list(counts = getCounts(sce)),
        colData = SummarizedExperiment::colData(sce),
        rowData = S4Vectors::DataFrame(
            simPICsanitizePeakAnnotationForRowData(
                as.data.frame(SummarizedExperiment::rowData(sce))
            )
        )
    )
    sce <- simPICfilterZeroLibraryCells(sce)
    sample.map <- simPICreadMicrogliaSampleMap(sample.map)

    if (!all(c("Sample", "Batch") %in% colnames(SummarizedExperiment::colData(sce)))) {
        stop("Microglia SCE must contain 'Sample' and 'Batch' columns.")
    }

    sce <- sce[, sce$Sample %in% sample.map$Sample, drop = FALSE]
    if (ncol(sce) == 0) {
        stop("No cells remain after subsetting to mapped Microglia samples.")
    }

    aligned <- simPICsubsetRowsByPeakAnnotation(sce, peak.annot)
    sce <- aligned$sce
    peak.annot <- aligned$peak.annot

    if (!("sample_batch" %in% colnames(SummarizedExperiment::colData(sce)))) {
        sce$sample_batch <- paste(sce$Sample, sce$Batch, sep = "_")
    }

    unit.sizes <- table(sce$sample_batch)
    keep.units <- names(unit.sizes)[unit.sizes > min.cells]

    if (!is.null(force.keep.batch)) {
        force.keep <- as.character(sce$Batch) == force.keep.batch
        if (!is.null(force.keep.pathology) &&
                "Pathology" %in% colnames(SummarizedExperiment::colData(sce))) {
            force.keep <- force.keep & as.character(sce$Pathology) == force.keep.pathology
        }
        keep.units <- union(
            keep.units,
            unique(as.character(sce$sample_batch[force.keep]))
        )
    }

    if (length(keep.units) == 0) {
        stop("No sample_batch units exceed min.cells in the Microglia example.")
    }

    sce <- sce[, sce$sample_batch %in% keep.units, drop = FALSE]
    unit.sizes <- sort(table(sce$sample_batch), decreasing = TRUE)
    big.unit <- names(unit.sizes)[1]
    bigcounts <- sce[, sce$sample_batch == big.unit, drop = FALSE]

    aggregated <- simPICaggregateMeanByGroup(sce, sce$sample_batch)

    keep.samples <- unique(as.character(SummarizedExperiment::colData(sce)$Sample))
    aligned.vcf <- simPICalignVcfToSamples(vcf, sample.map, keep.samples)
    peak.annot <- simPICharmonizePeakSeqnames(peak.annot, aligned.vcf$vcf)

    return(list(
        sce = sce,
        bigcounts = bigcounts,
        aggregated = aggregated,
        sample.map = aligned.vcf$sample.map,
        peak.annot = peak.annot[rownames(sce), , drop = FALSE],
        vcf = aligned.vcf$vcf
    ))
}
