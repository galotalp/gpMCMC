# (c) Copyright William J. Welch 2000-2006.
# 2000.11.14
# 2001.05.09: Test for ncol(x) > 1 in cor.fn.
# 2001.05.09: abs() inserted in h initialization for x[,1].
# 2009.03.13: as.matrix(x)
# 2010.08.10: cor.par coerced to matrix from vector if necessary

# Implemented for several explanatory variables - x is a matrix.

#' Correlation Matrix for Gaussian Process based on Power-Exponential Correlation function
#' (c) Copyright William J. Welch 2000-2006.
#'
#' @param x matrix/vector of explanatory variables
#' @param cor.par matrix of correlation parameters
#'
#' @return returns the correlation matrix based on power-exponential correlation structure
#' @export
#'
#' @examples n <- 5
#' x1 <- seq(-5,10,length.out = n)
#' x2 <- seq(0,15,length.out = n)
#'
#' data1 <- expand.grid(x1,x2)
#' x <- data1
#' # create hyperparameter matrix of thetas and alphas, alphas set to 0 indicated guassian correlation
#' d2 <- c(0.01,0.2,0,0)
#' cor.par <- data.frame(matrix(data = d2,nrow = dim(x)[2],ncol = 2))
#' names(cor.par) <- c("Theta.y","Alpha.y")
#'
#'
#' R <- cor.matrix(data1,cor.par) # obtain covariance matrix
#'
cor.matrix <- function(x, cor.par)
{


     # Make sure x is a matrix (not a dataframe).
     x <- as.matrix(x)


     # Make sure cor.par is a matrix (not a vector, i.e. in 1-d).
     if (is.vector(cor.par))
         cor.par <- matrix(cor.par, nrow = 1, ncol = length(cor.par))

     R <- matrix(nrow = nrow(x), ncol = nrow(x))


     #check to make sure dimensions are appropriate
     if(dim(x)[2]  != dim(cor.par)[1]) stop("mismatched dimensions")

     for (j in 1:ncol(R))
          R[, j] <- cor.fn(x, x[j, ], cor.par)

     return(R)
}

cor.fn <- function(x, x.vec, cor.par)
# Compute a vector of correlations between the rows of x
# and the vector x.vec.
# Assumes power-exponential correlation function.
# First column of cor.par contains theta's;
# second contains alpha's.
{
     theta <- cor.par[, 1]
     alpha <- cor.par[, 2]

     h <- theta[1] * abs(x[,1] - x.vec[1])^(2 - alpha[1])
     if (ncol(x) > 1)
          for (j in 2:ncol(x))
               h <- h + theta[j] * abs(x[,j] - x.vec[j])^(2 - alpha[j])

     return(exp(-h))
}




