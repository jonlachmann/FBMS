#######################################################
#
# Example 10 (Section 6.2): Zambia data set from the cAIC4 package
#
# Linear Mixed Model with Fractional Polynomials 
#
# Marginal Likelihood computed with lme4, INLA and with RTMB
#
# This is the valid version for the JSS Paper
#
#######################################################

library(tictoc)

library(devtools)
devtools::install_github("jonlachmann/GMJMCMC@FBMS", force=T, build_vignettes=F)

library(FBMS)
use.fbms = FALSE  



#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
#options(repos=c( inlabruorg = "https://inlabru-org.r-universe.dev", INLA = "https://inla.r-inla-download.org/R/testing", CRAN = "https://cran.rstudio.com") )
#install.packages("fmesher") 

library(lme4)
library(RTMB)
library(INLA)

#install.packages("cAIC4") 
#library(cAIC4)

data(Zambia, package = "cAIC4")

df <- as.data.frame(sapply(Zambia[1:5],scale))


transforms <- c("p0","p2","p3","p05","pm05","pm1","pm2","p0p0","p0p05","p0p1","p0p2","p0p3","p0p05","p0pm05","p0pm1","p0pm2")
probs <- gen.probs.gmjmcmc(transforms)
probs$gen <- c(1,1,0,1) # Modifications and interactions!

params <- gen.params.gmjmcmc(df)
params$feat$D <- 1   # Set depth of features to 1 (still allows for interactions)
params$loglik$r = 1/dim(df)[1]
params$feat$pop.max = 10

#specify indices for a random effect
params$loglik$dr = droplevels(Zambia$dr) # district ids for repeated measurements 


#estimator function with lme4

mixed.model.loglik.lme4 <- function (y, x, model, complex, params) 
{
  
  if (sum(model) > 1) {
    x.model = x[,model]
    data <- data.frame(y, x = x.model[,-1], dr = params$dr)

#    mm <- NULL
#    #importance with error handling for unstable libraries that one does not trust 100%
#    tryCatch({
#      mm <- lmer(as.formula(paste0("y ~ 1 +",paste0(names(data)[2:(dim(data)[2]-1)],collapse = "+"), "+ (1 | dr)")), data = data, REML = FALSE)
#    }, error = function(e) {
#      # Handle the error by setting result to NULL
#      mm <- NULL
#      # One can also print a message or log the error if needed
#      cat("An error in Estimation of MLIK occurred:", conditionMessage(e), "\n")
#    })
    mm <- lmer(as.formula(paste0("y ~ 1 +",paste0(names(data)[2:(dim(data)[2]-1)],collapse = "+"), "+ (1 | dr)")), data = data, REML = FALSE)
  } else{   #model without fixed effects
    data <- data.frame(y, dr = params$dr)
    mm <- lmer(as.formula(paste0("y ~ 1 + (1 | dr)")), data = data, REML = FALSE)
  }
    
    
   # logarithm of model prior
  if (length(params$r) == 0)  params$r <- 1/dim(x)[1]  # default value or parameter r
  lp <- log.prior(params, complex)
  
  mloglik <- as.numeric(logLik(mm))  -  0.5*log(length(y)) * (dim(data)[2] - 2) #Laplace approximation for beta prior
  
  return(list(crit = mloglik + lp, coefs = fixef(mm)))
}


#estimator function with INLA

params$loglik$INLA.num.threads = 10 # Number of threads used by INLA
#params$feat$keep.min = 0.2


mixed.model.loglik.inla <- function (y, x, model, complex, params) 
{
  if(sum(model)>1)
  {
    data1 = data.frame(y, as.matrix(x[,model]), params$dr)
    formula1 = as.formula(paste0(names(data1)[1],"~",paste0(names(data1)[3:(dim(data1)[2]-1)],collapse = "+"),"+ f(params.dr,model = \"iid\")"))
  } else
  {
    data1 = data.frame(y, params$dr)
    formula1 = as.formula(paste0(names(data1)[1],"~","1 + f(params.dr,model = \"iid\")"))
  }
  
  #to make sure inla is not stuck
  inla.setOption(inla.timeout=30)
  inla.setOption(num.threads=params$INLA.num.threads) 
  
  mod<-NULL
  #importance with error handling for unstable libraries that one does not trust 100%
  tryCatch({
    mod <- inla(family = "gaussian",silent = 1L,safe = F, data = data1,formula = formula1)
  }, error = function(e) {
    
    # Handle the error by setting result to NULL
    mod <- NULL
    
    # You can also print a message or log the error if needed
    cat("An error occurred:", conditionMessage(e), "\n")
  })
  
  # logarithm of model prior
  if (length(params$r) == 0)  params$r <- 1/dim(x)[1]  # default value or parameter r
  lp <- log.prior(params, complex)
  
  if(length(mod)<3||length(mod$mlik[1])==0) {
    return(list(crit = -10000 + lp,coefs = rep(0,dim(data1)[2]-2)))
  } else {
    mloglik <- mod$mlik[1]
    return(list(crit = mloglik + lp, coefs = mod$summary.fixed$mode))
  }
}


#estimator function with RTMB

params$loglik$nr_dr =  sum((table(Zambia$dr))>0)   #number of districts (that is number of different random intercepts)

mixed.model.loglik.rtmb <- function (y, x, model, complex, params) 
{
  
  z = model.matrix(y~params$dr) #Design matrix for random effect

  msize = sum(model)
  #Set up and estimate model
  dat = list(y = y, xm = x[,model], z = z)
  par = list(logsd_eps = 0,
             logsd_dr = 0,
             beta = rep(0,msize),
             u = rep(0,params$nr_dr))
  
  
  nll = function(par){
    getAll(par,dat)
    sd_eps = exp(logsd_eps)
    sd_dr = exp(logsd_dr)
    
    nll = 0
    #-log likelihood random effect
    nll = nll -  sum(dnorm(u, 0, sd_dr, log = TRUE))
    mu = as.vector(as.matrix(xm)%*%beta) + z%*%u
    nll <- nll - sum(dnorm(y, mu, sd_eps, log = TRUE))
    
#    ADREPORT(sd_dr) 
#    ADREPORT(sd_eps) 
    
    return(nll)
  }
  
  obj <- MakeADFun(nll , par, random = "u", silent = T )
#  obj <- MakeADFun(nll , par, random = "u")
  opt <- nlminb ( obj$par , obj$fn , obj$gr, control = list(iter.max = 10))

  # logarithm of model prior
  if (length(params$r) == 0)  params$r <- 1/dim(x)[1]  # default value or parameter r
  lp <- log.prior(params, complex)
  
#  if(length(beta)==0) {
#    return(list(crit = -10000 + lp,coefs = rep(0,dim(data1)[2]-2)))
#  } else {
  mloglik <- -opt$objective - 0.5*log(dim(x)[1])*msize
  return(list(crit = mloglik + lp, coefs = opt$par[-(1:2)]))
 # }
}



######################
#
# Compare runtime
#

set.seed(03052024)

#result <- gmjmcmc(data = df, loglik.pi = mixed.model.loglik.inla, transforms = transforms, probs = probs, params = params, P = 3)

#params$feat$pop.max = 10

tic()
if (use.fbms) {
  result1a <- fbms(data = df, family = "custom", loglik.pi = mixed.model.loglik.lme4, method = "gmjmcmc", 
                  transforms = transforms, N.init = 30, 
                  probs = probs, params = params, P=3)
} else { 
  result1a <- gmjmcmc(data = df, loglik.pi = mixed.model.loglik.lme4, 
                  transforms = transforms, N.init = 30, 
                  probs = probs, params = params, P = 3)
}
time.lme4 = toc()

plot(result1a)
summary(result1a, labels = names(df)[-1])


tic()
if (use.fbms) {
  result1b <- fbms(data = df, family = "custom", loglik.pi = mixed.model.loglik.inla, method = "gmjmcmc", 
                   transforms = transforms, N.init = 30, 
                   probs = probs, params = params, P=3)
} else { 
  result1b <- gmjmcmc(data = df, loglik.pi = mixed.model.loglik.inla, 
                    transforms = transforms, N.init = 30, 
                    probs = probs, params = params, P = 3)
}
time.inla = toc()

tic()
if (use.fbms) {
  result1c <- fbms(data = df, family = "custom", loglik.pi = mixed.model.loglik.rtmb, method = "gmjmcmc", 
                   transforms = transforms, N.init = 30, 
                   probs = probs, params = params, P=3)
} else { 
  result1c <- gmjmcmc(data = df, loglik.pi = mixed.model.loglik.rtmb, 
                      transforms = transforms, N.init = 30, 
                      probs = probs, params = params, P = 3)
}
time.rtmb = toc()
plot(result1c)
summary(result1c, labels = names(df)[-1])

c(time.lme4$callback_msg, time.inla$callback_msg, time.rtmb$callback_msg)

######################
#
# Analysis with lme4
#
#

set.seed(20062024)
params$feat$pop.max = 10
result2a <- gmjmcmc.parallel(runs = 40, cores = 10, data = df, loglik.pi = mixed.model.loglik.lme4, transforms = transforms, N.init=100, probs = probs, params = params, P = 25)

summary(result2a,tol = 0.05,labels=names(df)[-1])   


set.seed(21062024)
result2b <- gmjmcmc.parallel(runs = 120, cores = 40, data = df, loglik.pi = mixed.model.loglik.lme4, transforms = transforms, N.init=100, probs = probs, params = params, P = 25)

summary(result2b, labels = names(df)[-1])

summary(result2b, labels = names(df)[-1], pop = "all")
summary(result2b, labels = names(df)[-1], pop = "last")

plot(result2b)


set.seed(03072024)

result2c <- gmjmcmc.parallel(runs = 200, cores = 40, data = df, loglik.pi = mixed.model.loglik.lme4, transforms = transforms, N.init=100, probs = probs, params = params, P = 25)

summary(result2c, labels = names(df)[-1])
summary(result2c, labels = names(df)[-1], pop = "last")
summary(result2c, labels = names(df)[-1], pop = "all")
summary(result2c, labels = names(df)[-1], pop = "best")


summary(result2a, labels = names(df)[-1])
summary(result2b, labels = names(df)[-1])
summary(result2c, labels = names(df)[-1])


######################
#
# Analysis with INLA (Not used for manuscript, very long runtime)
#
#

set.seed(22052024)

params$loglik$INLA.num.threads = 1 # Number of threads used by INLA set to 1
result2a <- gmjmcmc.parallel(runs = 20, cores = 20, data = df, loglik.pi = mixed.model.loglik.inla, transforms = transforms, N.init=30, probs = probs, params = params, P = 25)

plot(result2a)
summary(result2a, labels = names(df)[-1])

#save.image("Ex9_Results2_parallel.RData")
#load("Ex9_Results_parallel.RData")

params$feat$check.col = F

set.seed(20062024)
params$loglik$INLA.num.threads = 1 # Number of threads used by INLA set to 1
result2b <- gmjmcmc.parallel(runs = 100, cores = 20, data = df, loglik.pi = mixed.model.loglik.inla, transforms = transforms, N.init=16, probs = probs, params = params, P = 15)

summary(result2b, labels = names(df)[-1])

