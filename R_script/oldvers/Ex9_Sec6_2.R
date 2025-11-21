#######################################################
#
# Example 9 (Section 6.2): Zambia data set from the cAIC4 package
#
# Linear Mixed Model with Fractional Polynomials, using only fbms
#
# Marginal Likelihood computed with lme4, INLA and with RTMB
#
# This is the valid version for the JSS Paper
#
#######################################################

library(tictoc)
library(FBMS)




#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
#options(repos=c( inlabruorg = "https://inlabru-org.r-universe.dev", INLA = "https://inla.r-inla-download.org/R/testing", CRAN = "https://cran.rstudio.com") )
#install.packages("fmesher") 

library(lme4)
library(RTMB)
library(INLA)

#install.packages("cAIC4") 
library(cAIC4)

data(Zambia, package = "cAIC4")
df <- as.data.frame(sapply(Zambia[1:5],scale))


transforms <- c("p0","p2","p3","p05","pm05","pm1","pm2","p0p0","p0p05","p0p1","p0p2","p0p3","p0p05","p0pm05","p0pm1","p0pm2")
probs <- gen.probs.gmjmcmc(transforms)
probs$gen <- c(1,1,0,1) # Modifications and interactions!

params <- gen.params.gmjmcmc(ncol(df) - 1)
params$feat$D <- 1   # Set depth of features to 1 (still allows for interactions)
params$feat$pop.max = 10


# function to estimate log posterior with lme4

mixed.model.loglik.lme4 <- function (y, x, model, complex, mlpost_params) 
{
  
  # logarithm of marginal likelihood (Laplace approximation)
  if (sum(model) > 1) {
    x.model = x[,model]
    data <- data.frame(y, x = x.model[,-1], dr = mlpost_params$dr)
    
    mm <- lmer(as.formula(paste0("y ~ 1 +",paste0(names(data)[2:(dim(data)[2]-1)],collapse = "+"), "+ (1 | dr)")), data = data, REML = FALSE)
  } else{   #model without fixed effects
    data <- data.frame(y, dr = mlpost_params$dr)
    mm <- lmer(as.formula(paste0("y ~ 1 + (1 | dr)")), data = data, REML = FALSE)
  }
  
  mloglik <- as.numeric(logLik(mm))  -  0.5*log(length(y)) * (dim(data)[2] - 2) #Laplace approximation for beta prior
  
  # logarithm of model prior
  if (length(mlpost_params$r) == 0)  mlpost_params$r <- 1/dim(x)[1]  # default value or parameter r
  lp <- log_prior(mlpost_params, complex)
  
  
  return(list(crit = mloglik + lp, coefs = fixef(mm)))
}


# function to estimate log posterior with INLA

mixed.model.loglik.inla <- function (y, x, model, complex, mlpost_params) 
{
  if(sum(model)>1)
  {
    data1 = data.frame(y, as.matrix(x[,model]), mlpost_params$dr)
    formula1 = as.formula(paste0(names(data1)[1],"~",paste0(names(data1)[3:(dim(data1)[2]-1)],collapse = "+"),"+ f(mlpost_params.dr,model = \"iid\")"))
  } else
  {
    data1 = data.frame(y, mlpost_params$dr)
    formula1 = as.formula(paste0(names(data1)[1],"~","1 + f(mlpost_params.dr,model = \"iid\")"))
  }
  
  #to make sure inla is not stuck
  inla.setOption(inla.timeout=30)
  inla.setOption(num.threads=mlpost_params$INLA.num.threads) 
  
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
  if (length(mlpost_params$r) == 0)  mlpost_params$r <- 1/dim(x)[1]  # default value or parameter r
  lp <- log_prior(mlpost_params, complex)
  
  if(length(mod)<3||length(mod$mlik[1])==0) {
    return(list(crit = -10000 + lp,coefs = rep(0,dim(data1)[2]-2)))
  } else {
    mloglik <- mod$mlik[1]
    return(list(crit = mloglik + lp, coefs = mod$summary.fixed$mode))
  }
}


# function to estimate log posterior with RTMB

mixed.model.loglik.rtmb <- function (y, x, model, complex, mlpost_params) 
{
  z = model.matrix(y~mlpost_params$dr) #Design matrix for random effect
  
  msize = sum(model)
  #Set up and estimate model
  dat = list(y = y, xm = x[,model], z = z)
  par = list(logsd_eps = 0,
             logsd_dr = 0,
             beta = rep(0,msize),
             u = rep(0,mlpost_params$nr_dr))
  
  nll = function(par){
    getAll(par,dat)
    sd_eps = exp(logsd_eps)
    sd_dr = exp(logsd_dr)
    
    nll = 0
    #-log likelihood random effect
    nll = nll -  sum(dnorm(u, 0, sd_dr, log = TRUE))
    mu = as.vector(as.matrix(xm)%*%beta) + z%*%u
    nll <- nll - sum(dnorm(y, mu, sd_eps, log = TRUE))
    
    return(nll)
  }
  obj <- MakeADFun(nll , par, random = "u", silent = T )
  opt <- nlminb ( obj$par , obj$fn , obj$gr, control = list(iter.max = 10))
  
  # logarithm of model prior
  if (length(mlpost_params$r) == 0)  mlpost_params$r <- 1/dim(x)[1]  # default value or parameter r
  lp <- log_prior(mlpost_params, complex)
  
  mloglik <- -opt$objective - 0.5*log(dim(x)[1])*msize
  return(list(crit = mloglik + lp, coefs = opt$par[-(1:2)]))
}



######################
#
# Compare runtime
#

set.seed(03052024)

tic()
result1a <- fbms(formula = z ~ 1+., data = df, transforms = transforms,
                 method = "gmjmcmc",probs = probs, params = params, P=3, N = 30,
                 family = "custom", loglik.pi = mixed.model.loglik.lme4,
                 model_prior = list(r = 1/dim(df)[1]), 
                 extra_params = list(dr = droplevels(Zambia$dr)))
time.lme4 = toc()


tic()
result1b <- fbms(formula = z ~ 1+., data = df, transforms = transforms,
                 method = "gmjmcmc",probs = probs, params = params, P=3, N = 30,
                 family = "custom", loglik.pi = mixed.model.loglik.inla,
                 model_prior = list(r = 1/dim(df)[1]), 
                 extra_params = list(dr = droplevels(Zambia$dr), 
                                     INLA.num.threads = 10))
time.inla = toc()

tic()
result1c <- fbms(formula = z ~ 1+., data = df, transforms = transforms,
                 method = "gmjmcmc",probs = probs, params = params, P=3, N = 30,
                 family = "custom", loglik.pi = mixed.model.loglik.rtmb,
                 model_prior = list(r = 1/dim(df)[1]), 
                 extra_params = list(dr = droplevels(Zambia$dr), 
                                    nr_dr =  sum((table(Zambia$dr))>0)))
time.rtmb = toc()
plot(result1c)
summary(result1c, labels = names(df)[-1])

c(time.lme4$callback_msg, time.inla$callback_msg, time.rtmb$callback_msg)

######################
#
# Analysis with lme4
#
#


# Only this one is actually reported in the manuscript
set.seed(20062024)
params$feat$pop.max = 10

result2a <- fbms(formula = z ~ 1+., data = df, transforms = transforms,
                 probs = probs, params = params, P=25, N = 100,
                 method = "gmjmcmc.parallel", runs = 40, cores = 40,
                 family = "custom", loglik.pi = mixed.model.loglik.lme4,
                 model_prior = list(r = 1/dim(df)[1]), 
                 extra_params = list(dr = droplevels(Zambia$dr)))
  
summary(result2a,tol = 0.05,labels=names(df)[-1])   


set.seed(21062024)

result2b <- fbms(formula = z ~ 1+., data = df, transforms = transforms,
                 probs = probs, params = params, P=25, N = 100,
                 method = "gmjmcmc.parallel", runs = 120, cores = 40,
                 family = "custom", loglik.pi = mixed.model.loglik.lme4,
                 model_prior = list(r = 1/dim(df)[1]), 
                 extra_params = list(dr = droplevels(Zambia$dr)))

summary(result2b, labels = names(df)[-1])

summary(result2b, labels = names(df)[-1], pop = "all")
summary(result2b, labels = names(df)[-1], pop = "last")

plot(result2b)


set.seed(03072024)

result2c <- fbms(formula = z ~ 1+., data = df, transforms = transforms,
                 probs = probs, params = params, P=25, N = 100,
                 method = "gmjmcmc.parallel", runs = 200, cores = 40,
                 family = "custom", loglik.pi = mixed.model.loglik.lme4,
                 model_prior = list(r = 1/dim(df)[1]), 
                 extra_params = list(dr = droplevels(Zambia$dr)))

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

# Number of threads used by INLA set to 1
result2aI <- fbms(formula = z ~ 1+., data = df, transforms = transforms,
                  probs = probs, params = params, P=25, N = 100,
                  method = "gmjmcmc.parallel", runs = 40, cores = 40,
                  family = "custom", loglik.pi = mixed.model.loglik.inla,
                  model_prior = list(r = 1/dim(df)[1]), 
                  extra_params = list(dr = droplevels(Zambia$dr), INLA.num.threads = 1))

plot(result2aI)
summary(result2aI, labels = names(df)[-1])



params$feat$check.col = F

set.seed(20062024)
# Number of threads used by INLA set to 1
result2bI <- fbms(formula = z ~ 1+., data = df, transforms = transforms,
                  probs = probs, params = params, P=25, N = 100,
                  method = "gmjmcmc.parallel", runs = 120, cores = 40,
                  family = "custom", loglik.pi = mixed.model.loglik.inla,
                  model_prior = list(r = 1/dim(df)[1]), 
                  extra_params = list(dr = droplevels(Zambia$dr), INLA.num.threads = 1))

plot(result2bI)
summary(result2bI, labels = names(df)[-1])