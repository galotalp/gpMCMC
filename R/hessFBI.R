# this version of HESS FBI includes the derivative w.r.t the exponential prior
# this function takes non-transformed theta parameters and takes the derivative with respect to log-transformed theta parameters.

	# you need to do the following before running this function:
	# Read in the following matrices:
	#cor.par <- read.mtx("corpar.mtx") these are theta values, not tau values
	#x <- read.mtx("x.mtx")
	#y <- read.mtx("y.mtx")
	#f <- read.mtx("F.mtx")
	# Remove "Term" column
	# cor.par <- cor.par[, -1]

#' Hessian Matrix for Hyper-parameters in a Guassian Process Posterior distribution with an Exponential Prior
#'
#' this version of HESS FBI includes the derivative w.r.t the exponential prior
#' this function takes non-transformed theta parameters and takes the derivative with respect to log-transformed theta parameters.
#'
#' @param x covariate matrix/ vector x
#' @param y response vector
#' @param f regression model, must be a matrix. if constant should be a vector (as.matrix) of 1s of length n, where n is the number of data points in x
#' @param cor.par matrix of theta and alpha parameters for power-exponential model, includes two columns, the first for theta parameters and the second for alpha parameters.  For Guassian correlation structure, alpha parameters can be initialized as 0.
#' @param lambda parameter for exponential prior, defaults to 0.1
#'
#' @return returns a matrix of partial second derivatives w.r.t tau = log(theta) correlation parameters.  Taking the negative inverse of this hessian matrix provides the fisher approximation of the covariance
#' @export
#'
#' @examples
#' n <- 5
#' x1 <- seq(-5,10,length.out = n)
#' x2 <- seq(0,15,length.out = n)
#' x <- expand.grid(x1,x2)
#' d2 <- c(0.01,0.2,0,0)
#' cor.par <- data.frame(matrix(data = d2,nrow = dim(x)[2],ncol = 2))
#' names(cor.par) <- c("Theta.y","Alpha.y")
#'
#' R <- cor.matrix(x,cor.par) # obtain covariance matrix
#' L <- chol(R)
#' z <- as.vector(rnorm(n^2))
#' y <- t(L)%*%z
#'
#' f <- as.matrix(rep(1,n^2))
#' hess <- hessFBI(x,y,f,cor.par)
hessFBI <- function(x,y,f,cor.par, lambda = 0.1){

	f <- as.matrix(f)
	y <- as.matrix(y)


	R <- cor.matrix(x, cor.par)

	U <- chol(R)
	UT <- t(U)

	t <- 0
	n <- dim(x)[1]
	k <- dim(x)[2]
	m <- dim(f)[2]
	nd <- n - m

	out <- mat.or.vec(k,k)
	yTilde <- forwardsolve(UT,y)
	yDoubleT <- backsolve(U,yTilde)
	FTilde <- forwardsolve(UT,f)
	q <- qr(FTilde,tol = 1e-07 ,useLAPACK = TRUE)
	Q1 <- qr.Q(q, complete = FALSE)
	R1 <- qr.R(q, complete = FALSE)

	Q1T <- t(Q1)

	yb <- Q1T%*%yTilde

	varHat <- (t(yTilde)%*%yTilde - t(yb)%*%yb)/nd

	f1 <- Q1%*%yb
	f2 <- yTilde - f1
	f3 <- backsolve(U,f2)

	yf1 <- Q1T%*%y
	yf2 <- backsolve(R1,yf1)
	yf3 <- f%*%yf2
	yfb <- y - yf3
	yfb2 <- forwardsolve(UT,yfb)
	yfb3 <- backsolve(U, yfb2)


	dr_dt <- mat.or.vec(n,n)
	d2r_dt <- mat.or.vec(n,n)
	dr_dt2 <- mat.or.vec(n,n)
	dv_dt <- mat.or.vec(k,1)


	V1 <- mat.or.vec(n,n)
	V2 <- mat.or.vec(n,n)
	V3 <- mat.or.vec(n,n)
	V4 <- mat.or.vec(n,n)
	V5 <- mat.or.vec(n,n)
	V6 <- mat.or.vec(n,n)
	V7 <- mat.or.vec(n,n)
	Term1 <- mat.or.vec(m,n)
	Term2 <- mat.or.vec(n,n)
	Term3 <- mat.or.vec(n,m)
	Cfin <- mat.or.vec(m,m)

	sum <- 0.0
	r1 <- mat.or.vec(n,n)
	r2 <- mat.or.vec(n,n)
	r3 <- mat.or.vec(n,n)
	r4 <- mat.or.vec(n,n)
	r5 <- mat.or.vec(n,n)
	r6 <- mat.or.vec(n,n)
	ya <- mat.or.vec(n,1)
	y1 <- mat.or.vec(n,1)
	y2 <- mat.or.vec(n,1)
	I <- diag(n)
	qq <- Q1%*%Q1T

	sigma1 <- 0.0
	sigma2 <- 0.0
	dsigma <- 0.0
	trace <- 0.0
	trace2 <- 0.0

	for(t in 1:k){
		for(i in 1:n){
			for(j in 1:n){
				if(i == j){dr_dt[i,j] <- 0.0}
				else {dr_dt[i,j] <- R[i,j]*cor.par[t,1]*(-1.0)*(abs(x[i,t]-x[j,t]))^(2.0 - cor.par[t,2])}
				}
			}
		i <- 0.0
		j <- 0.0

		ya <- dr_dt%*%f3
		dv_dt[t] <- t(f3)%*%ya
		dv_dt[t] <- -dv_dt[t]/nd

		}

	for(i in 1:k){
		sigma1 <- dv_dt[i]

		for(t in 1:n){
			for(l in 1:n){
				if(t==l){dr_dt[t,l] <- 0.0}

				else{dr_dt[t,l] <- R[t,l]*(-1.0)*cor.par[i,1]*(abs(x[t,i]-x[l,i]))^(2.0 - cor.par[i,2])}
			}
		}

		for(j in 1:k){
			sigma2 <- dv_dt[j]
			for(t in 1:n){
				for(l in 1:n){
					if(t==l){
						dr_dt2[t,l] <- 0.0
						d2r_dt[t,l] <- 0.0
						}

					else{
						dr_dt2[t,l] <- R[t,l]*cor.par[j,1]*(-1.0)*(abs(x[t,j]-x[l,j]))^(2.0 - cor.par[j,2])
						if(i==j){
							d2r_dt[t,l] <- R[t,l]*cor.par[i,1]*(abs(x[t,i]-x[l,i]))^(2.0 - cor.par[i,2])*(cor.par[i,1]*(abs(x[t,i]-x[l,i]))^(2.0 - cor.par[i,2]) - 1)
							}
						else{
							d2r_dt[t,l] <- R[t,l]*cor.par[i,1]*cor.par[j,1]*(abs(x[t,i]-x[l,i]))^(2.0 - cor.par[i,2])*(abs(x[t,j]-x[l,j]))^(2.0 - cor.par[j,2])
							}

						}
					}
				}

			r1 <- forwardsolve(UT,dr_dt)
			r2 <- forwardsolve(UT,dr_dt2)
			r1 <- t(r1)
			r3 <- I - qq
			r3 <- r1%*%r3%*%r2
			r3 <- 2*r3 - d2r_dt
			y1 <- r3%*%yfb3
			dsigma <- t(yfb3)%*%y1

			r4 <- forwardsolve(UT,dr_dt2)
			r4 <- backsolve(U,r4)
			r4 <- dr_dt%*%r4
			r4 <- forwardsolve(UT,r4)
			r4 <- backsolve(U,r4)

			r5 <- forwardsolve(UT,d2r_dt)
			r5 <- backsolve(U,r5)
			r6 <- r5 - r4
			trace <- 0.0

			for(t in 1:n) {trace = trace + r6[t,t]}


			Term3 <- backsolve(U,FTilde)
			Term1 <- backsolve(U, Q1)
			Term1 <- backsolve(R1, t(Term1))

			V5 <- dr_dt2
			V6 <- dr_dt

			V4 <- d2r_dt

			V3 <- forwardsolve(UT,dr_dt)
			V3 <- backsolve(U, V3)
			V3 <- V5%*%V3

			V2 <- forwardsolve(UT,dr_dt2)
			V2 <- backsolve(U,V2)
			V2 <- V6%*%V2

			V7 <- V2 - V3

			V1 <- forwardsolve(UT,dr_dt)
			V1 <- Q1T%*%V1
			V1 <- Q1%*%V1
			V1 <- backsolve(U,V1)

			V1 <- V5%*%V1

			Term2 <- V1 + V2 - V3 + V4

			Cfin <- Term1%*%Term2%*%Term3
			trace2 <- 0.0
			for(t in 1:m) {trace2 = trace2 + Cfin[t,t]}


			out[i,j] = ((0.5*nd)/(varHat^2))*dv_dt[i]*dv_dt[j] - dsigma/(2.0*varHat) -0.5*trace + 0.5*trace2
			if(i==j){
				out[i,j] = out[i,j] -lambda*(cor.par[i,1])
			}

			}
	}
    return(out)
}
