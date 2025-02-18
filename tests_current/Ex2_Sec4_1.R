#######################################################
#
# Example 2 (Section 4.1):
#
# Simulated data without any nonlinearities
#
# This is the valid version for the JSS Paper
#
#######################################################

# Logical to decide whether to perform analysis with fbms function
# If FALSE then gmjmcmc or gmjmcmc.parallel function is used
use.fbms <-  FALSE  
stronger.singal <- FALSE

library(mvtnorm)
library(FBMS)



n <- 100  # sample size
p <- 20   # number of covariates
p.vec <- 1:p


k <- 5    #size of the data generating model

set.seed(1002)

correct.model <- sample(p.vec, k)
beta.k <- 1*ifelse(stronger.singal,10,1) + rnorm(k)/2   # Coefficents of the correct submodel

beta <- c(rep(0, p))
beta[correct.model] <- beta.k

x <- rmvnorm(n, rep(0, p))
y <- x %*% beta    + rnorm(n)
X <- as.matrix(x)

y<-scale(y)
X<-scale(X)/sqrt(n)


df <- as.data.frame(cbind(y, X))

correct.model
beta.k


to3 <- function(x) x^3
transforms <- c("sigmoid","sin_deg","exp_dbl","p0","troot","to3")

set.seed(123)
if (use.fbms) {
  result <- result <- fbms(data = df, method = "gmjmcmc", 
                           transforms = transforms, P = 40)
} else {
  result <- gmjmcmc(df, gaussian.loglik, gaussian.loglik.alpha, transforms, P = 40)
}
summary(result)


set.seed(123)
if (use.fbms) {
  result2 <- result <- fbms(data = df, method = "gmjmcmc", transforms = transforms, 
                           N.init = 1000, N.final = 5000, P = 40)
} else {
  result2 <- gmjmcmc(df, transforms = transforms, 
                     N.init = 1000, N.final = 5000, P = 40)
}
summary(result2)




########################################################
#
#  Model which includes no non-linear effects
#
#


set.seed(123)

probs.lin <- gen.probs.mjmcmc()

params.lin <- gen.params.mjmcmc(df)
params.lin$loglik$r <- 1/dim(df)[1] 

#to set variance to unknown uncomment below
#params.lin$loglik$var <- "unknown"

if (use.fbms) {
  result.lin <- fbms(data = df, N = 5000)
} else {
  result.lin <- mjmcmc(df, N = 5000, probs = probs.lin, 
                       params =  params.lin)
  
}


plot(result.lin)
summary(result.lin)

correct.model
beta.k


#The default value of  N = 100  does not lead to the same result here
set.seed(123)

if (use.fbms) {
  result.lindef <- fbms(data = df)
} else {
  result.lindef <- mjmcmc(df)
}


plot(result.lindef)
summary(result.lindef)


###############################################################################