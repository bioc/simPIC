#' Reads counts data from an RDS file and performs simPIC simulation and 
#' estimation.
#'
#' @details This code is an example demonstrating how to use the provided test
#' data to get the default parameters as in \code{\link{simPICcount}}
#' 
#' @seealso \code{\link{newsimPICcount}}, \code{\link{simPICestimate}}
#' 
#' @name testdata
#' @rdname testdata

options(max.print = 30)

counts <- readRDS("inst/extdata/test.rds")
new <- newsimPICcount()
new <- simPICestimate(counts, pm.distr ="weibull")
print(new)
