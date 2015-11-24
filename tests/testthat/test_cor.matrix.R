test_that("cor.matrix generates correlation matrix of appropriate dimension",{
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

  expect_true(dim(R)[1] == dim(data1)[1])
  expect_true(dim(R)[2] == dim(data1)[1])
  expect_true(is.numeric(R))

})


test_that("cor.matrix expects error at incorrect dimensions of cor.par",{
  n <- 10
  x1 <- seq(-5,10,length.out = n)
  x2 <- seq(0,15,length.out = n)

  data1 <- expand.grid(x1,x2)
  x <- data1
  # create hyperparameter matrix of thetas and alphas, alphas set to 0 indicated guassian correlation
  d2 <- c(0.01,0.2,0.4,0,0,0)
  cor.par <- data.frame(matrix(data = d2,nrow = 3,ncol = 2))
  names(cor.par) <- c("Theta.y","Alpha.y")


  expect_error(cor.matrix(data1,cor.par), "mismatched dimensions")# obtain covariance matrix


})
