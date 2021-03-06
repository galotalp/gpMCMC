% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fitGauss.R
\name{fitGauss}
\alias{fitGauss}
\title{Runs gasp code to do gaussian process inference using covariates x and y and a constant regression model and a Gaussian correlation function.  Returns the fitted theta parameters}
\usage{
fitGauss(x, y)
}
\arguments{
\item{x}{covariates}

\item{y}{response}
}
\value{
fitted theta parameters
}
\description{
Runs gasp code to do gaussian process inference using covariates x and y and a constant regression model and a Gaussian correlation function.  Returns the fitted theta parameters
}
\examples{
n <- 10
x1 <- seq(-5,10,length.out = n)
x2 <- seq(0,15,length.out = n)
x <- expand.grid(x1,x2)
x <- as.matrix(x)
d2 <- c(0.01,0.2,0,0) #here we set the theta parameters to be 0.01 and 0.2.
cor.par <- data.frame(matrix(data = d2,nrow = dim(x)[2],ncol = 2))
names(cor.par) <- c("Theta.y","Alpha.y")

R <- cor.matrix(x,cor.par) # obtain covariance matrix
L <- chol(R)
z <- as.matrix(rnorm(n^2))
y <- L\%*\%z

fitGauss(x,y)
}

