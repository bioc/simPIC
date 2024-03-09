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
#' This function converts a sparse matrix into a SingleCellExperiment(SCE)
#' object.
#'
#' @param sparse_data A sparse matrix containing count data, where rows are
#' peaks and columns represent cells.
#' @return A SingleCellExperiment(SCE) object with the sparse matrix stored
#' in the "counts" assay.
#' @importFrom SingleCellExperiment SingleCellExperiment
#' @importFrom Matrix Matrix
convert_to_SCE <- function(sparse_data) {
    sce <- SingleCellExperiment(assays = list(counts = sparse_data))
    rownames(sce) <- paste0("Peak", seq_len(nrow(sparse_data)))
    colnames(sce) <- paste0("Cell", seq_len(ncol(sparse_data)))
    return(sce)
}
