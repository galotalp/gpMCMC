#### MCMC 3, multivariate normal distributed steps, each step is one p-dimensional vector from a multivariate normal distribution
# code for sampling gaussian correlation parameters in a gaussian process

#read.mtx function for reading in matrix files
# source("readmtx.R")
# library(MASS)

#here we read in a matrix of the maximum likelihood estimators for the correlation coefficients, theta and alpha.
# maximum likelihood estimators for theta values are used as initial estimates for the metropolis hastings algorithm

# cor.par <- read.mtx("corpar.mtx")
# cor.par <- cor.par[,-1]
# cor.par2 <- cor.par

bceMCMC_mvrnorm <-function(nmcmc,burn,thin,x,y,reg,step, priortheta){

  #here we read in a matrix of the maximum likelihood estimators for the correlation coefficients, theta and alpha.
  # maximum likelihood estimators for theta values are used as initial estimates for the metropolis hastings algorithm

  cor.par <- read.mtx("corpar.mtx")
  cor.par <- cor.par[,-1]
  cor.par2 <- cor.par

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
          phi.cond<-log(phi)+mvrnorm(n = 1,rep(0, p),diag(p)*phi.w)  #currently using normal distribution with mean 0 and variance 1
          if(all(exp(phi.cond) > 0)){
            d <- 1
          }
        }

        phi.cond <- exp(phi.cond)
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

