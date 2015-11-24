test_that("gpMCMC generates numeric values",{
  nsamp <- 100
  burn <- 200
  thin <- 10
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

  gp <- gpMCMC(nsamp,burn,thin,x,y,step = 1)
  m1 <-  mean(gp$samples[,2]) #these means should be similar to the theta parameters set above
  m2 <- mean(gp$samples[,1])
  expect_is(m1,"numeric")
  expect_is(m2,"numeric")

})

test_that("gpMCMC generates a sample matrix of specific dimensions",{
  nsamp <- 100
  burn <- 200
  thin <- 10
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

  gp <- gpMCMC(nsamp,burn,thin,x,y,step = 1)

  expect_true(dim(gp$samples)[1] == nsamp)

})



test_that("gpMCMC generates an acceptance vector of specific dimensions",{
  nsamp <- 100
  burn <- 200
  thin <- 10
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

  gp <- gpMCMC(nsamp,burn,thin,x,y,step = 1,method = "uniform")

  expect_true(length(gp$acceptance) == 4)

})


test_that("gpMCMC generates an acceptance vector of specific dimensions, again for normally distributed steps",{
  nsamp <- 100
  burn <- 200
  thin <- 10
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

  gp <- gpMCMC(nsamp,burn,thin,x,y,step = 1,method = "normal")

  expect_true(length(gp$acceptance) == 4)

})





