#### MCMC 2, normally distributed Metropolis algorithm proposal step lengths.
# code for sampling gaussian correlation parameters in a gaussian process
##### MCMC code for sampling gaussian correlation parameters in a gaussian process

#' Univariate, Normally distributed step length MCMC function, meant only to be used by gpMCMC function
#'
#' @param nmcmc number of MCMC samples to be generated before thinning and burning
#' @param burn number of mcmc samples to burn
#' @param thin keep one of every 'thin' samples
#' @param x covariates
#' @param y response
#' @param reg only option currently is "constant"
#' @param step step length for mcmc
#' @param priortheta only option currently is "Exp"
#'
#' @return returns a list containing mcmc.ma (samples) and accept (acceptance rates)
#' @export
#'
#' @examples
#'
#' nsamp <- 100
#' burn <- 200
#' thin <- 10
#'
#' n <- 10
#' x1 <- seq(-5,10,length.out = n)
#' x2 <- seq(0,15,length.out = n)
#' x <- expand.grid(x1,x2)
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
#' gp <- bceMCMC_norm(1000,10,10,x,y,reg = "constant",step =0.1, priortheta = "Exp")
#' mean(gp$mcmc.ma[,2]) #these means should be similar to the theta parameters set above
#' mean(gp$mcmc.ma[,1])
bceMCMC_norm <-function(nmcmc,burn,thin,x,y,reg,step, priortheta){

  ddd <- c(rep(1,dim(x)[2]),rep(0,dim(x)[2]))
  cor.par <- data.frame(matrix(data = ddd,nrow = dim(x)[2],ncol = 2))
  names(cor.par) <- c("Theta.y","Alpha.y")
  cor.par2 <- cor.par

  p<-ncol(x)
  j=0

  #final mcmc trials
  mcmc.ma<-matrix(nrow=(nmcmc-burn)/thin,ncol=ncol(x))

  #we use the following vector to count acceptances
  accept<-matrix(c(rep(0,p+2)),nrow=1,ncol=p+2,byrow=T)

  #initial guesses
  phi<-c(rep(0.1,p))
  #phi<-cor.par[,1]

  #step length
  phi.w<-c(rep(step,p))




  if(reg=="constant"){

    for(i in 1:nmcmc){

      for(k in 1:p){

        phi.cond<-phi
        d <- 0
        while(d == 0){
          phi.cond[k]<-log(phi[k])+(rnorm(1))*phi.w[k]  #currently using normal distribution with mean 0 and variance 1
          if(exp(phi.cond[k])>0){
            d <- 1
          }
        }
	phi.cond[k] <- exp(phi.cond[k])

        if(phi.cond[k]>0){
          phi_cond=phi.cond
          phi_or=phi
          #com.phi<-log.post1.constant(x,y,phi_cond, priortheta)$logpost-log.post1.constant(x,y,phi_or, priortheta)$logpost
          cor.par[,1] <- phi_cond
          cor.par2[,1] <- phi_or
          com.phi <- log_posterior(x,y,as.matrix(rep(1,dim(x)[1])),cor.par, prior = "Exp") - log_posterior(x,y,as.matrix(rep(1,dim(x)[1])),cor.par2, prior = "Exp")
          u<-runif(1)
          if(log(u)<com.phi){
            phi<-phi.cond
            accept[1,(2+k)]=accept[1,(2+k)]+1
          }
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


