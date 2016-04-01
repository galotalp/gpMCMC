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

library(MASS)

fitGauss(x,y)



d <- gpMCMC(100,100,10,x,y)
mean(d$samples[,1])
mean(d$samples[,2])

