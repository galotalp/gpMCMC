RandomNumberSeed = 500

# Input files:
X < x1.mtx
Y < y1.mtx
RegressionModel < constmod.mtx

CorrelationFamily = PowerExponential
RandomError = No
Alpha.Max = 0
# MLE
CriticalLogLikelihoodDifference = 0
LogLikelihoodTolerance = 0.001
Tries = 3
FitCriterion = posterior
Fit

CrossValidate

# Output files:
StochasticProcessModel > corpar2.mtx


Quit
