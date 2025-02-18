#######################################################
#
# Example 12 (Section 6.4):
#
# Subsampling
#
# Heart Disease Health Indicators Dataset”
#
# This is the valid version for the JSS Paper
#
#######################################################

#install.packages("tictoc")
library(tictoc)

library(devtools)
devtools::install_github("jonlachmann/GMJMCMC@FBMS", force=T, build_vignettes=F)
#install.packages("FBMS")
library(FBMS)
#library(devtools)
#devtools::install_github("jonlachmann/irls.sgd", force=T, build_vignettes=F)
library(irls.sgd)



setwd("/home/florian/FBMS/")

df = read.csv2(file = "heart_disease_health_indicators_BRFSS2015.csv",sep = ",",dec = ".")

summary(df)
dim(df)

#number of observations in the data

n = dim(df)[1] 

#number of covariates

p = dim(df)[2] - 1   

params <- gen.params.gmjmcmc(data = df)
params$loglik$r = 0.5 
params$loglik$subs = 0.01


transforms <- c("sigmoid","pm1","p0","p05","p2","p3")
probs <- gen.probs.gmjmcmc(transforms)

logistic.posterior.bic.irlssgd <- function (y, x, model, complex, params) 
{
  if (!is.null(params$crit)) {
    mod <- glm.sgd(x[,model], y, binomial(), sgd.ctrl = list(start=params$coefs, subs=params$subs, maxit=10, alpha=0.00008, decay=0.99, histfreq=10))
    mod$deviance <- get_deviance(mod$coefficients, x[,model], y, binomial())
    mod$rank <- length(mod$coefficients)
  } else {
    mod <- irls.sgd(as.matrix(x[,model]), y, binomial(),
                  irls.control=list(subs=params$subs, maxit=20, tol=1e-7, cooling = c(1,0.9,0.75), expl = c(3,1.5,1)),
                  sgd.control=list(subs=params$subs, maxit=250, alpha=0.001, decay=0.99, histfreq=10))
  }
  
  # logarithm of marginal likelihood
  mloglik <- -mod$deviance / 2 - 0.5 * log(length(y)) * (mod$rank - 1)
    
  # logarithm of model prior
  if (length(params$r) == 0)  params$r <- 1/dim(x)[1]  # default value or parameter r
  lp <- log.prior(params, complex)
  crit <- mloglik + lp

  if (!is.null(params$crit) && params$crit > crit) {
    return(list(crit = params$crit, coefs = params$coefs))
  }
  
  return(list(crit = crit, coefs = mod$coefficients))
}



############################
#
# Testing runtime
#
############################

set.seed(100001)
tic()
# subsampling analysis
tmp1 <- gmjmcmc(df, logistic.posterior.bic.irlssgd, transforms = transforms, 
                params = params, P = 2, sub  = T)
time1 = toc() 

set.seed(100002)
tic()
# regular analysis
tmp2 <- gmjmcmc(df, logistic.loglik, transforms = transforms, 
                params = params, P = 2)
time2 = toc() 

c(time1, time2)

############################
#
# More serious analysis
#
############################


# with subsampling

set.seed(100003)

tic()
result <- gmjmcmc.parallel(runs = 10,cores = 10, data = df,
                           loglik.pi = logistic.posterior.bic.irlssgd, 
                           transforms = transforms, params = params,  P = 3, sub = T)
time3 = toc() 

summary(result)

# without subsampling


set.seed(100004)

tic()
result1a <- gmjmcmc.parallel(runs = 10,cores = 10, data = df,
                           loglik.pi = logistic.loglik, 
                           transforms = transforms, params = params,  P = 3)
time4 = toc() 

summary(result1a)



############################
#
# Probably too ambitious analysis
#
############################

# with subsampling

set.seed(100005)
tic()
result2 <- gmjmcmc.parallel(runs = 40,cores = 40, data = df,
                           loglik.pi = logistic.posterior.bic.irlssgd, 
                           transforms = transforms, params = params,  P = 10, sub = T)
time5 = toc() 
summary(result2)



# regular analysis


set.seed(100006)

tic()
result2a <- gmjmcmc.parallel(runs = 40,cores = 40, data = df,
                             loglik.pi = logistic.loglik, 
                             transforms = transforms, params = params,  P = 10)
time6 = toc() 


summary(result2a)



############################################################################

C = cor(df, use = "everything",
        method = "spearman")

corrplot::corrplot(C)

apply((abs(C - diag(diag(C)))), 2, max)

save.image("Ex11_Results.RData")

