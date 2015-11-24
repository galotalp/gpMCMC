#' For writing out matrix files to be read by gasp software
#'
#' @param x matrix to be written
#' @param filename intended name of matrix file, should end with ".mtx"
#' @param row.name optional vector indicating row-names, by default set to row-number
#' @param col.name optional vector indicating column-names, by default set to x1, x2, etc
#' @param page.width default = 80 characters
#'
#' @return returns nothing
#' @export
#'
#' @examples
#'
#' x <- c(1:3)
#' write.mtx(x,"x.mtx")
write.mtx <- function(x, filename, row.name = dimnames(x)[[1]],
   col.name = dimnames(x)[[2]], page.width = 80)
# Write matrix / data frame x to filename in blocked ACED/GaSP format.
# Row and column labels are preserved if present.
# If there are no row labels, then 1, 2, ... will be used.
# If there are no column labels, then x1, x2, ... will be
# used, unless another name is given for col.name.
# 1997.01.24: replaced with new version.
# 2001.11.08: Renamed write.mx (from write.mat).
# 2002.10.17: Renamed write.mtx (from write.mx).
# 2007.03.04: Handling of row and column names changed. Argument row.name
#             introduced (set to all elements to "" for no row labels);
#             argumment col.name now allows a vector.
{
     # Make sure x is a matrix.
     x <- as.matrix(x)

     nrows    <- nrow(x)
     ncolumns <- ncol(x)

     # Row labels
     if (!is.null(row.name))
          row.labels <- row.name
     else
          row.labels <- dimnames(x)[[1]]
     if (is.null(row.labels))
          row.labels <- c(1:nrows)
     if (sum(row.labels != "") == 0)
          row.labels <- c("----", "", "----", row.labels)
     else
          row.labels <- c("----", "Case", "----", row.labels)

     # Column labels
     if (!is.null(col.name))
          col.labels <- col.name
     else
          col.labels <- dimnames(x)[[2]]
     if (is.null(col.labels))
          col.labels <- paste("x", 1:ncolumns, sep="")

     spaces <- c(rep("----", ncolumns))

     mat <- rbind(spaces, col.labels, spaces, x)
     dimnames(mat) <- list(row.labels, c(rep(" ", ncol(mat))))

     sink(filename)
     options(length=nrow(mat)+10)
     # Prevent repeats of column labels.
     temp <- options(length = nrow(mat) + 1, width = page.width)

     print(mat,quote=F)
     sink()
     cat(filename, ": ", nrow(mat) - 3, " rows and ", ncol(mat),
               " columns written.\n", sep = "")

     # Restore length.
     options(temp)
}

