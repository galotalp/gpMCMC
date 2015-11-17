#### MCMC 4, multivariate normal distributed steps, each step is one p-dimensional vector from a multivariate normal distribution,
# in this case the multivariate normal distribution is the laplace approximation to the variance
# code for sampling gaussian correlation parameters in a gaussian process

#read.mtx function for reading in matrix files
# source("readmtx.R")
# source("hessFBI.R")
# library(MASS)

#here we read in a matrix of the maximum likelihood estimators for the correlation coefficients, theta and alpha.
# maximum likelihood estimators for theta values are used as initial estimates for the metropolis hastings algorithm
#
# cor.par <- read.mtx("corpar.mtx")
# cor.par <- cor.par[,-1]
# cor.par2 <- cor.par
#f <- read.mtx("F.txt")




bceMCMC_mvrnorm_2 <-function(nmcmc,burn,thin,x,y,reg,step, priortheta){

  #here we read in a matrix of the maximum likelihood estimators for the correlation coefficients, theta and alpha.
  # maximum likelihood estimators for theta values are used as initial estimates for the metropolis hastings algorithm

  cor.par <- read.mtx("corpar.mtx")
  cor.par <- cor.par[,-1]
  cor.par2 <- cor.par
  f <- as.matrix(rep(1,dim(x)[1]))


  hess <- hessFBI(x,y,f,cor.par)
 cov1 <- -solve(hess)

  p<-ncol(x)
  j=0



  #final mcmc trials
  mcmc.ma<-matrix(nrow=(nmcmc-burn)/thin,ncol=ncol(x))

  #we use the following vector to count acceptances
  accept<- 0

  #initial guesses
  phi<-c(rep(0.1,p))
  #phi<-cor.par[,1]

  #step length
  phi.w<-c(rep(step^2,p)/p)



  if(reg=="constant"){

    for(i in 1:nmcmc){

      phi.cond<-phi
      d <- 0
      while(d == 0){
        phi.cond <- log(phi)+mvrnorm(n = 1,rep(0, p),cov1)
        if(all(exp(phi.cond) > 0)){
          d <- 1
        }
      }
      phi.cond <- exp(phi.cond)
#
# this part in between here is just for testing
#     while(d == 0){
#          phi.cond<- phi+mvrnorm(n = 1,rep(0, p),cov1)
#           if(all(phi.cond > 0)){
#             d <- 1
#            }
#          }
# this part above here is just for debugging

      if(all(phi.cond > 0)){
        phi_or=phi
        phi_cond=phi.cond
        cor.par[,1] <- phi_cond
        cor.par2[,1] <- phi_or
        com.phi <- log.posterior(x,y,as.matrix(rep(1,dim(x)[1])),cor.par, prior = "Exp") - log.posterior(x,y,as.matrix(rep(1,dim(x)[1])),cor.par2, prior = "Exp")
        u<-runif(1)
        if(log(u)<com.phi){
          phi<-phi.cond
          accept <- accept + 1
        }
      }

      if(i>burn&&((i-burn)%%thin==0)){
        j=j+1
        mcmc.ma[j,]=phi
      }



      #if(i>burn&&((i-burn)%%thin==0)){
      #j=j+1
      #res[j,]=pred1.constant(x,y,xtest1,mcmc.ma[i,3:(p+2)], priortheta, priorsigma)$res
      #v.term2[j,]=pred1.constant(x,y,xtest1,mcmc.ma[i,3:(p+2)], priortheta, priorsigma)$v.term2
      #}

      if ((i%%(0.1*nmcmc))==0){
        print(c(paste("reg=", reg), paste("priortheta=", priortheta) ,i/nmcmc))
        #print(c('prior=1',i/nmcmc))
      }

    }
    #mcmc.ma<-mcmc.ma[,-(1:2)]
    m<-list(mcmc.ma=mcmc.ma, accept = accept)
    return(m)
  }

}

