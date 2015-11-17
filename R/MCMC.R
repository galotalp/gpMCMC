##### MCMC code for sampling gaussian correlation parameters in a gaussian process

#read.mtx function for reading in matrix files
# source("readmtx.R")

#here we read in a matrix of the maximum likelihood estimators for the correlation coefficients, theta and alpha.
# maximum likelihood estimators for theta values are used as initial estimates for the metropolis hastings algorithm

# cor.par <- read.mtx("corpar.mtx")
# cor.par <- cor.par[,-1]
# cor.par2 <- cor.par

bceMCMC<-function(nmcmc,burn,thin,x,y,reg,step, priortheta){

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
			  phi.cond[k]<-log(phi[k])+(runif(1)-0.5)*phi.w[k]
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
				com.phi <- log.posterior(x,y,as.matrix(rep(1,dim(x)[1])),cor.par, prior = "Exp") - log.posterior(x,y,as.matrix(rep(1,dim(x)[1])),cor.par2, prior = "Exp")
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



