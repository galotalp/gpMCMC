test_that("log_posterior generates a value for the posterior likelihood",{
  n <- 5
  x1 <- seq(-5,10,length.out = n)
  x2 <- seq(0,15,length.out = n)

  data1 <- expand.grid(x1,x2)
  x <- data1
  # create hyperparameter matrix of thetas and alphas, alphas set to 0 indicated guassian correlation
  d2 <- c(0.01,0.2,0,0)
  cor.par <- data.frame(matrix(data = d2,nrow = dim(x)[2],ncol = 2))
  names(cor.par) <- c("Theta.y","Alpha.y")


  R <- cor.matrix(data1,cor.par) # obtain covariance matrix
  L <- chol(R)
  z <- as.vector(rnorm(n^2))
  y <- t(L)%*%z

  logpost <- log_posterior(data1,y,as.matrix(rep(1,n^2)),cor.par,prior = "Exp")
  expect_is(logpost,"numeric")

})

test_that("log_posterior is greater for theta values far from the posterior mode",{
  n <- 5
  x1 <- seq(-5,10,length.out = n)
  x2 <- seq(0,15,length.out = n)

  data1 <- expand.grid(x1,x2)
  x <- data1
  # create hyperparameter matrix of thetas and alphas, alphas set to 0 indicated guassian correlation
  d2 <- c(0.01,0.2,0,0)
  cor.par <- data.frame(matrix(data = d2,nrow = dim(x)[2],ncol = 2))
  names(cor.par) <- c("Theta.y","Alpha.y")


  R <- cor.matrix(data1,cor.par) # obtain covariance matrix
  L <- chol(R)
  z <- as.vector(rnorm(n^2))
  y <- t(L)%*%z


  d3 <- c(rnorm(1,4,2),rnorm(1,4,2),0,0)
  cor.par2 <- data.frame(matrix(data = d3,nrow = dim(x)[2],ncol = 2))

  logpost <- log_posterior(data1,y,as.matrix(rep(1,n^2)),cor.par,prior = "Exp")
  logpost2 <- log_posterior(data1,y,as.matrix(rep(1,n^2)),cor.par2,prior = "Exp")
  expect_gt(logpost,logpost2)
})


test_that("Only Accepts Matrices of appropriate dimensions",{
  n <- 5
  x1 <- seq(-5,10,length.out = n)
  x2 <- seq(0,15,length.out = n)

  data1 <- expand.grid(x1,x2)
  x <- data1
  # create hyperparameter matrix of thetas and alphas, alphas set to 0 indicated guassian correlation
  d2 <- c(0.01,0.2,0,0)
  cor.par <- data.frame(matrix(data = d2,nrow = dim(x)[2],ncol = 2))
  names(cor.par) <- c("Theta.y","Alpha.y")


  R <- cor.matrix(data1,cor.par) # obtain covariance matrix
  L <- chol(R)
  z <- as.vector(rnorm(n^2))
  y <- t(L)%*%z

  expect_error(log_posterior(c(14,13),y,c(1,1,1,1,1,1),cor.par,prior = "Exp"))
})


