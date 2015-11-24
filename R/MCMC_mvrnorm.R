
# code for sampling gaussian correlation parameters in a gaussian process
#' This is an MCMC function with a multivariate normal proposal density that is only meant to be used by the gpMCMC master function
#'
#' #### MCMC 3, multivariate normal distributed steps, each step is one p-dimensional vector from a multivariate normal distribution
#'
#' @param nmcmc Number of mcmc samples to generate, before burning and thinning
#' @param burn number of samples to burn
#' @param thin keep only one of every 'thin' samples
#' @param x predictors
#' @param y response
#' @param reg currently only option is "constant"
#' @param step $(step^2)/p$ is multiplied by the identity matrix, which is then used as the covariance for a multivariate normal proposal density in the Metropolis Algorithm
#' @param priortheta only currently only option is "Exp", "Higs" and "none" will also be implemented
#'
#' @return returns a list containing mcmc.ma (samples) and accept (acceptance rates)
#' @export
#'
#' @examples
#' nsamp <- 100
#' burn <- 200
#' thin <- 10
#'
#' n <- 10
#' x1 <- seq(-5,10,length.out = n)
#' x2 <- seq(0,15,length.out = n)
#' x <- expand.grid(x1,x2)
#' x <- as.matrix(x)
#' d2 <- c(0.01,0.2,0,0) #here we set the theta parameters to be 0.01 and 0.2.
#' # These are the modes of the distribution that we will sample from using MCMC
#' cor.par <- data.frame(matrix(data = d2,nrow = dim(x)[2],ncol = 2))
#' names(cor.par) <- c("Theta.y","Alpha.y")
#'
#' R <- cor.matrix(x,cor.par) # obtain covariance matrix
#' L <- chol(R)
#' z <- as.matrix(rnorm(n^2))
#' y <- L%*%z
#'
#' gp <- bceMCMC_mvrnorm(1000,10,10,x,y,reg = "constant",step =1, priortheta = "Exp")
#' mean(gp$mcmc.ma[,2]) #these means should be similar to the theta parameters set above
#' mean(gp$mcmc.ma[,1])
bceMCMC_mvrnorm <-function(nmcmc,burn,thin,x,y,reg,step, priortheta){
  ddd <- c(rep(1,dim(x)[2]),rep(0,dim(x)[2]))
  cor.par <- data.frame(matrix(data = ddd,nrow = dim(x)[2],ncol = 2))
  names(cor.par) <- c("Theta.y","Alpha.y")
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
          phi.cond<-log(phi)+MASS::mvrnorm(n = 1,rep(0, p),diag(p)*phi.w)  #currently using normal distribution with mean 0 and variance 1
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
          com.phi <- log_posterior(x,y,as.matrix(rep(1,dim(x)[1])),cor.par, prior = "Exp") - log_posterior(x,y,as.matrix(rep(1,dim(x)[1])),cor.par2, prior = "Exp")
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

