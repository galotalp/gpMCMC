





#' master MCMC function
#'
#' @param nsamp number of MCMC samples to be generated
#' @param burn number of samples to be burned, must be less than nsamp
#' @param thin only include one of every 'thin' samples.  If thin is 10, only include one in every 10 samples
#' @param x covariate matrix/ vector x
#' @param y response vector
#' @param reg regression model, defaults to "constant"
#' @param step step length for MCMC, defaults to 1
#' @param priortheta prior, defaults to Exponential prior
#' @param method proposal distribution, defaults to "uniform"
#' @param plot defaults to FALSE
#'
#' @return returns a list containing the mcmc samples (mcmc_samples) and the acceptance rates for each covariate (acceptance)
#' @export
#'
#' @examples
gpMCMC <-function(nsamp,burn,thin,x,y,reg = "constant",step = 1, priortheta = "Exp", method = "uniform", plot = FALSE){

  if(thin <= 0){
    nmcmc <- nsamp + burn
  }
  else{
    nmcmc <- nsamp*thin + burn
  }

  mcmc_samples <- data.frame(matrix(data = NA, nrow = nsamp, ncol = dim(x)[2]))
  acceptance <- 0

  if(reg == "constant"){
    if(priortheta == "Exp" || priortheta == "exp"){
       if(method == "uniform"){
         re <- bceMCMC(nmcmc,burn,thin,x,y,reg="constant",step = step, priortheta="Exp")
         mcmc_samples <- re$mcmc.ma
         acceptance <- re$accept
       }
       else if(method == "normal"){
         re <- bceMCMC_norm(nmcmc,burn,thin,x,y,reg="constant",step = step, priortheta="Exp")
         mcmc_samples <- re$mcmc.ma
         acceptance <- re$accept
       }
       else if(method == "mvrnorm"){
         re <- bceMCMC_mvrnorm(nmcmc,burn,thin,x,y,reg="constant",step = step, priortheta="Exp")
         mcmc_samples <- re$mcmc.ma
         acceptance <- re$accept
       }
       else if(method == "mvrnorm2"){
         re<-bceMCMC_mvrnorm_2(nmcmc,burn,thin,x,y,reg="constant",step = step, priortheta="Exp")
         mcmc_samples <- re$mcmc.ma
         acceptance <- re$accept
       }
    }
#     if(priortheta == "Gem"){
#       break
#     }
#     if(priortheta == "Jeffreys"){
#       break
#     }
  }
names(mcmc_samples) <- names(x)
if(plot == FALSE){
 m <- list(samples = mcmc_samples, acceptance = acceptance)
 return(m)
}
else{
  m <- list(samples = mcmc_samples, acceptance = acceptance)
  t1s <- mat.or.vec(200,dim(x)[2])


  for(j in 1:dim(x)[2]){
    t1s[,j] <-   seq(log(cor.par[j,1]) - 2*0.50*abs(log(cor.par[j,1])),log(cor.par[j,1]) + 2*0.495*abs(log(cor.par[j,1])),10^(-2)*abs(log(cor.par[j,1])))
  }


  cor.par2 <- cor.par
  corpar <- cor.par
  names1 <- colnames(x)


  #LL1 to LL6 are obviously specific to the g-protein dataset and this functionality won't be applicable for other problems
#   c <- contour(t1s[,2],t1s[,1],LL1,xlab = names1[2], ylab = names1[1],main = "log-posterior contour with mcmc samples")
#   points(log(mcmc_samples[,2]),log(mcmc_samples[,1]),pch = 20)
#   show(c)
#
#   c <- contour(t1s[,3],t1s[,1],LL2,xlab = names1[3], ylab = names1[1],main = "log-posterior contour with mcmc samples")
#   points(log(mcmc_samples[,3]),log(mcmc_samples[,1]),pch = 20)
#   show(c)
#
#   c <- contour(t1s[,4],t1s[,1],LL3,xlab = names1[4], ylab = names1[1],main = "log-posterior contour with mcmc samples")
#   points(log(mcmc_samples[,4]),log(mcmc_samples[,1]),pch = 20)
#   show(c)
#
#   c <- contour(t1s[,3],t1s[,2],LL4,xlab = names1[3], ylab = names1[2],main = "log-posterior contour with mcmc samples")
#   points(log(mcmc_samples[,3]),log(mcmc_samples[,2]),pch = 20)
#   show(c)
#
#   c <- contour(t1s[,4],t1s[,2],LL5,xlab = names1[4], ylab = names1[2],main = "log-posterior contour with mcmc samples")
#   points(log(mcmc_samples[,4]),log(mcmc_samples[,2]),pch = 20)
#   show(c)
#
#   c <- contour(t1s[,4],t1s[,3],LL6,xlab = names1[4], ylab = names1[3],main = "log-posterior contour with mcmc samples")
#   points(log(mcmc_samples[,4]),log(mcmc_samples[,3]),pch = 20)
#   show(c)
  return(m)
}
}








