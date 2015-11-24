# gpMCMC

Simulating from a Gaussian Process Posterior distribution using MCMC and Fast Bayesian Inference

This package is currently under developement.

 The package provides the means for sampling from Posterior Distributions of Gaussian Processes that use the Gaussian Correlation function. This package also provides a function for computing the Fisher approximation to the covariance of this posterior likelihood distribution. 
 
 
## To be added:

  1. Additional priors: currently the only prior being used is the Exponential Prior
  2. Ability to make calls to gasp software so that posterior inference can be done in C. This will make functions mcmc_mvrnorm_2 and fitGauss functional (these functions are currently not functional)
  3. A plotting mechanism for plotting the samples against pairs of covariates.
  4. An "FBI-Sampling" method, for sampling from the Fischer Approximated distribution of hyper-parameters

## Tests:

Currently tests are implemented for checking all working MCMC sampling function.  MCMC_mvrnorm_2 will currently not provide adequate results, and therefore the "mvrnorm2" option for methods in the master gpMCMC function will also provide false results.  

