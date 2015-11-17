read.mtx <- function(filename, sep = "")
# Read a matrix in ACED/GaSP format from filename and return a 
# data frame.
# This version can handle blocked matrices, the row and column 
# names are preserved, and a mixture of numeric and character
# columns is allowed.
# 1995.05.04: Now handles empty matrices.
# 1995.05.12: row.labels error corrected.
# 2001.11.08: Renamed read.mx (from read.mat) and returns a data frame.
# 2002.02.26: sep added as argument (passed to scan);  
#             use sep = "," for comma-delimited files.  
# 2002.02.20: argument matrix.mode removed; columns with all numerical
#             values (or NA's) are coerced to numeric.
# 2002.10.17: Renamed read.mtx (from read.mx).
# (c) Copyright William J. Welch 1995-2002.
{
     z <- scan(filename, what = character(), sep = sep)

     # Set up "---", "----", etc. and "___", "____", etc.
     max.len <- max(nchar(z))
     dash <- character(max.len - 2)
     underscore <- character(max.len - 2)
     for (i in 3 : max.len)
     {
          dash[i-2] <- paste(rep("-", times = i), collapse = "")
          underscore[i-2] <- paste(rep("_", times = i), collapse = "")
     }

     # Indices of z corresponding to dashes (or underscores).
     is.dash <- match(z, c(dash, underscore), nomatch = 0)
     dash.index <- (1:length(is.dash))[is.dash > 0]

     # Add a terminating dash.
     dash.index[length(dash.index) + 1] <- length(z) + 1
     
     block.num <- 0

     # While more blocks to read.
     while (length(dash.index) > 2 &&
               dash.index[3] - dash.index[2] > 1)
     {
          block <- matrix(z[(dash.index[2] + 1) : (dash.index[3] - 1)],
                    ncol = dash.index[2] - dash.index[1] - 1,
                    byrow = TRUE, dimnames = list(NULL,
                    z[(dash.index[1] + 1) : (dash.index[2] - 1)]))
          block.num <- block.num + 1

          col1.name <- dimnames(block)[[2]][1]
          if (col1.name == "case" || col1.name == "Case" ||
                    col1.name == "CASE")
          {
               row.labels <- as.character(block[, 1])
               block <- block[, -1, drop = F]
          }
          else
               row.labels <- NULL

          # Add this block to x.
          if (block.num == 1)
               x <- block
          else
               x <- cbind(x, block)

          dash.index <- dash.index[-c(1,2)]
     }

     if (block.num == 0 || length(dash.index) != 1)
     {
          cat(filename, ": Returning NULL.\n", sep = "")
          return(NULL)
     }
     else
     {
          dimnames(x) <- list(row.labels, dimnames(x)[[2]])
          x <- data.frame(x)

          # Coerce numerical columns to numeric mode.
          for (j in 1:ncol(x))
          {
               xj <- x[, j]

               # Suppress warning messages from coercing.
               options(warn = -1)
               xj.lev <- as.numeric(levels(xj))
               options(warn = 0)

               if (sum(is.na(xj.lev)) == 0)
                    # All levels are numeric.
                    x[, j] <- xj.lev[xj]
          }

          cat(filename, ": ", nrow(x), " rows and ", ncol(x),
               " columns read.\n", sep = "")

          return(x)
     }
}

