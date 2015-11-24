#' Runs gasp code to do gaussian process inference using covariates x and y and a constant regression model and a Gaussian correlation function.  Returns the fitted theta parameters
#'
#' @param x covariates
#' @param y response
#'
#' @return fitted theta parameters
#' @export
#'
#' @examples
#' n <- 10
#' x1 <- seq(-5,10,length.out = n)
#' x2 <- seq(0,15,length.out = n)
#' x <- expand.grid(x1,x2)
#' x <- as.matrix(x)
#' d2 <- c(0.01,0.2,0,0) #here we set the theta parameters to be 0.01 and 0.2.
#' cor.par <- data.frame(matrix(data = d2,nrow = dim(x)[2],ncol = 2))
#' names(cor.par) <- c("Theta.y","Alpha.y")
#'
#' R <- cor.matrix(x,cor.par) # obtain covariance matrix
#' L <- chol(R)
#' z <- as.matrix(rnorm(n^2))
#' y <- L%*%z
#'
#' fitGauss(x,y)
fitGauss <- function(x,y){

#   write.mtx(x,paste(system.file("bin", package = "proj2"),"/x1.mtx", sep = ""))
#   write.mtx(y,paste(system.file("bin", package = "proj2"),"/y1.mtx", sep = ""))
#   system(paste(system.file("bin", "gasp", package = "proj2"),system.file("bin", "fit.gsp", package = "proj2")))

  corp <- read.mtx(system.file("bin", "corpar2.mtx", package = "proj2"))

  thetas <- corp[,2]
  return(thetas)

}