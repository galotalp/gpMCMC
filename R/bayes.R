# necessary to execute function
# source("cor.R")
# source("readmtx.R")

# x <- read.mtx("x.mtx")
# y <- read.mtx("y.mtx")
# f <- read.mtx("F.mtx")
# cor.par <- read.mtx("corpar.mtx")

# f <- as.matrix(f)
# y <- as.matrix(y)

# cor.par <- cor.par[, -1]

#takes non-transformed theta parameters

log.posterior <- function(x, y, f, cor.par, prior)
{
     R <- cor.matrix(x, cor.par)


     n <- nrow(R)
     k <- dim(f)[2]



    U1<-try(chol(R),silent=T)
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


