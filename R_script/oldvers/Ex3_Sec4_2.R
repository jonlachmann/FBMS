#######################################################
#
# Example 3 (Section 4.2):
#
# Simulated data with interactions, using only fbms
#
# This is the valid version for the JSS Paper
#
#######################################################

library(mvtnorm)
library(FBMS)

n <- 100  # sample size
p <- 20   # number of covariates

# Model:  
# X1: Pure Main effect
# X2 : X3: Pure interaction effect
# X4 * X5: Main effects plus interaction effect


set.seed(1003)

x = rmvnorm(n, rep(0, p))
X <- as.matrix(x)
X <- scale(X)/sqrt(n)

y <- (1.2 * x[,1] + 1.5 * x[,2]* x[,3] - x[,4] + 1.1*x[,5] - 1.3 * x[,4]*x[,5])+ rnorm(n)
y<-scale(y)

df <- data.frame(y = y, X)


transforms <- c("")
probs <- gen.probs.gmjmcmc(transforms)
probs$gen <- c(1,0,0,1)            #Include interactions and mutations

####################################################
#
# single thread analysis (two different runs)
#
####################################################

set.seed(123)
result <- fbms(data = df, method = "gmjmcmc", transforms = transforms, 
                 probs = probs)
summary(result)


set.seed(123)
result2 <- fbms(data = df, method = "gmjmcmc", transforms = transforms, 
                 N = 1000, probs = probs, P=40)
summary(result2, tol = 0.01)


####################################################
#
# multiple thread analysis
#
####################################################


set.seed(123)

  result_parallel <- fbms(data = df, method = "gmjmcmc.parallel", transforms = transforms,
                          runs = 40, cores = 40,
                          probs = probs, P=25)

summary(result_parallel, tol = 0.01)



# Using longer more iterations of MJMCMC chains does not change results substantially
set.seed(123)

result_parallel2 <- fbms(data = df, method = "gmjmcmc.parallel", transforms = transforms,
                           runs = 40, cores = 10, N=1000, N.final=2000,
                           probs = probs, P=25)
summary(result_parallel2, tol = 0.01)

#summary(result_parallel2, pop = "all", tol = 0.01)


