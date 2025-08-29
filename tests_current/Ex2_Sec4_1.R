#######################################################
#
# Example 2 (Section 4.1):
#
# Simulated data without any nonlinearities, only using fbms function
#
# This is the valid version for the JSS Paper
#
#######################################################

#library(devtools)
#devtools::install_github("jonlachmann/FBMS@v1_arxiv", force=T, build_vignettes=F)

library(mvtnorm)
library(FBMS)



n <- 100  # sample size
p <- 20   # number of covariates
p.vec <- 1:p


k <- 5    #size of the data generating model


correct.model <- 1:k
beta.k <- (1:5)/5   # Coefficents of the correct submodel

beta <- c(rep(0, p))
beta[correct.model] <- beta.k

set.seed(123)

x <- rmvnorm(n, rep(0, p))
y <- x %*% beta    + rnorm(n)
X <- as.matrix(x)

y<-scale(y)
X<-scale(X)/sqrt(n)


df <- as.data.frame(cbind(y, X))
colnames(df) <- c("Y", paste0("X", seq_len(ncol(df) - 1)))

correct.model
beta.k

########################################################
#
#  Models with non-linear effects (gmjmcmc)
#
#

to3 <- function(x) x^3
transforms <- c("sigmoid","sin_deg","exp_dbl","p0","troot","to3")

set.seed(1)
  result <- fbms(data = df, method = "gmjmcmc", transforms = transforms)
  summary(result)
  plot(result)

   
set.seed(2)
  result2 <- fbms(data = df, method = "gmjmcmc", transforms = transforms, 
                            N = 1000, P = 40)
  summary(result2, tol = 0.1)
  plot(result)



########################################################
#
#  Model which includes no non-linear effects (mjmcmc)
#
#

  # The default value of  N = 1000 works relatively well here. 
  set.seed(1)
  result.lindef <- fbms(data = df)
  summary(result.lindef)
  plot(result.lindef)
  
  # Check that this is actually the default
  set.seed(1)
  result.lin <- fbms(data = df, N = 1000)
  summary(result.lin)
  plot(result.lin)
  



