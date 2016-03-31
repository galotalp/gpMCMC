library(gpMCMC)

n <- 10
x1 <- seq(-5,10,length.out = n)
x2 <- seq(0,15,length.out = n)
x <- expand.grid(x1,x2)
x <- as.matrix(x)
d2 <- c(2,1,0,0) #here we set the theta parameters to be 0.01 and 0.2.
cor.par <- data.frame(matrix(data = d2,nrow = dim(x)[2],ncol = 2))
names(cor.par) <- c("Theta.y","Alpha.y")
f <- as.matrix(rep(1,n^2))
R <- cor.matrix(x,cor.par) # obtain covariance matrix
L <- chol(R)
z <- as.matrix(rnorm(n^2))
y <- L%*%z


fitGauss(x,y)

library(MASS)
init <- rep(0.1, dim(x)[2])

logpost_wrapper <- function(val){
	d2 <- c(val[1],val[2],0,0) # FIX THIS: MUST generalize to all dimensions
	cor.par <- data.frame(matrix(data = d2,nrow = dim(x)[2],ncol = 2))
	names(cor.par) <- c("Theta.y","Alpha.y")
	return(log_posterior(x, y, f, cor.par, prior = "Exp"))
}

fit2 <- optim(init,logpost_wrapper, lower = -10, upper = 10)

d <- gpMCMC(100,100,10,x,y)
mean(d$samples[,1])
mean(d$samples[,2])

