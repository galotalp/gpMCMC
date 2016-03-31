#takes non-transformed theta parameters
#' log.posterior for gaussian process for x,y, and f, at the given correlation hyperparameters specified in cor.par
#'
#' @param x covariate matrix/vector x
#' @param y response vector y
#' @param f regression model, must be a matrix. If constant, should be a vector (as.matrix) of 1s of length n, where n is the number of data points in x
#' @param cor.par matrix of theta and alpha parameters for power-exponential model, includes two columns, the first for theta parameters and the second for alpha parameters.  For Guassian correlation structure, alpha parameters can be initialized as 0.
#' @param prior prior for log-posterior, options are "Exp" for exponential prior and "Hig" for Higdon's prior
#'
#' @return returns the log-posterior
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
#' L <- chol(R)
#' z <- as.vector(rnorm(n^2))
#' y <- t(L)%*%z
#'
#' logpost <- log_posterior(data1,y,as.matrix(rep(1,n^2)),cor.par,prior = "Exp")
log_posterior <- function(x, y, f, cor.par, prior)
{
     R <- cor.matrix(x, cor.par)


     n <- nrow(R)
     k <- dim(f)[2]

    U1<-try(chol(R),silent=T)  # catch error, if error is caught return a huge number.
	F1=mat.or.vec(n,1)+1
	S=backsolve(U1,F1,transpose=T)
	S1=crossprod(S)
	G=backsolve(U1,y,transpose=T)
	A=crossprod(S,G)
	##mu.hat
	beta.hat=solve(S1,A)

	B=y-beta.hat*F1
	B1=backsolve(U1,B,transpose=T)
	B2=crossprod(B1)
	#sigma_hat=(1/(n-k))t(y-Fbeta_hat)(R)^-1(y-Fbeta_hat)
	sigma.sq.hat=(1/(n-k))*B2








     U <- chol(R)
	 UT <- t(U)

     FTilde <- forwardsolve(UT,f)
     FTildeT <- t(FTilde)

     B <- FTildeT%*%FTilde


     # det expects a Matrix.  Remember "library Matrix".
     #R <- Matrix(as.vector(R), nrow = n, ncol = n)

     # Put in all constants - especially n.
     # Similarly gasp code.
     # with prior


     if(prior=="Hig"){
     	tem<-0.125*exp(-0.25*cor.par[,1])*(1-exp(-0.25*cor.par[,1]))^(-0.5)
     	logprior<-sum(log(tem))

    	log.lik <- logprior +  as.numeric(-0.5 * ((n-k) * log(sigma.sq.hat) + determinant(R, logarithm = T)$modulus + determinant(B, logarithm = T)$modulus) )
     }
     else if(prior=="Exp"){
     	lambda <- 0.1
     	tem <- (log(lambda) - lambda*cor.par[,1])
     	logprior<-sum(tem)

    	log.lik <- logprior +  as.numeric(-0.5 * ((n-k) * log(sigma.sq.hat) + determinant(R, logarithm = T)$modulus + determinant(B, logarithm = T)$modulus) )
     }
     else{
     	log.lik <- as.numeric(-0.5 * ((n-k) * log(sigma.sq.hat) + determinant(R, logarithm = T)$modulus + determinant(B, logarithm = T)$modulus) )

     }

     return(log.lik)
}


