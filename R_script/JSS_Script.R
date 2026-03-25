#################################################
#
# Example 1 (Section 3):
#
# Kepler Example with the most recent database update
#
# Basic introduction of the FBMS package
#
##################################################



#install.packages("devtools")
library(devtools)
install_github("jonlachmann/FBMS@jss_v2", force=T, build_vignettes=F)

library(FBMS)

data(exoplanet)

train.indx <- 1:500
df.train = exoplanet[train.indx, ]
df.test = exoplanet[-train.indx, ]


to3 <- function(x) x^3
transforms <- c("sigmoid","sin_deg","exp_dbl","p0","troot","to3")


####################################################
#
# single thread analysis (default values, Section 3.1)
#
####################################################


set.seed(123)

result.default <- fbms(formula = semimajoraxis ~ 1 + . , data = df.train, 
                       method = "gmjmcmc", transforms = transforms)


####################################################
#
# Choosing a different prior (more iterations, Section 3.3)
#
####################################################


set.seed(234)

result.BIC <- fbms(formula = semimajoraxis ~ 1 + . , data = df.train, 
                       method = "gmjmcmc", transforms = transforms,
                       beta_prior = list(type = "Jeffreys-BIC", Var = "unknown"))


set.seed(345)

result.EB <- fbms(formula = semimajoraxis ~ 1 + . , data = df.train, 
                      method = "gmjmcmc", transforms = transforms,
                      beta_prior = list(type = "EB-global", a = 1))

####################################################
#
# single thread analysis (more iterations, Section 3.4)
#
####################################################


set.seed(123)

result.P50 <- fbms(data = df.train, method = "gmjmcmc", transforms = transforms,
                     P = 50, N = 1000, N.final = 5000)

 
####################################################
#
# multiple thread analysis (Section 3.5)
#
####################################################

set.seed(123)

result.parallel <- fbms(data = df.train, method = "gmjmcmc.parallel", transforms = transforms,
                          runs = 40, cores = parallel::detectCores()-1, P = 25)


####################################################
#
# Inspection of Results (Section 3.6)
#
####################################################

######################
# summary

summary(result.default)
summary(result.default, pop = "all", labels = paste0("x",1:length(df.train[,-1])))


summary(result.P50)
summary(result.P50, pop = "best", labels = paste0("x",1:length(df.train[,-1])))
summary(result.P50, pop = "last", labels = paste0("x",1:length(df.train[,-1])))
summary(result.P50, pop = "last", tol = 0.01, labels = paste0("x",1:length(df.train[,-1])))
summary(result.P50, pop = "all")

summary(result.parallel)
library(tictoc)
tic()
summary(result.parallel, tol = 0.01, pop = "all",data = df.train)
toc()




######################
# plot

pdf("result.pdf") 
plot(result.default)
dev.off()

plot(result.default)



pdf("result.P50.pdf") 
plot(result.P50)
dev.off()

plot(result.P50)



pdf("result.parallel.pdf") 
plot(result.parallel)
dev.off()

plot(result.parallel)


####################################################
#
# Prediction (Section 3.7)
#
####################################################


#preds <- predict(result.default, df.test[,-1], link = function(x) x)
preds <-  predict(result.default, df.test[,-1])

str(aggr(preds))



rmse.default <- sqrt(mean((predmean(preds) - df.test$semimajoraxis)^2))

pdf("prediction.pdf") 
plot(predmean(preds), df.test$semimajoraxis)
dev.off()

plot(predmean(preds), df.test$semimajoraxis)






###############################


#preds.P50 = predict(result.P50, df.test[,-1], link = function(x) x)  
preds.P50 = predict(result.P50, df.test[,-1])  
rmse.P50 <-  sqrt(mean((predmean(preds.P50) - df.test$semimajoraxis)^2))

pdf("prediction.P50.pdf") 
plot(predmean(preds.P50), df.test$semimajoraxis)
dev.off()

plot(predmean(preds.P50), df.test$semimajoraxis)



###############################


preds.multi <- predict(result.parallel , df.test[,-1], link = function(x) x)
rmse.parallel <- sqrt(mean((predmean(preds.multi) - df.test$semimajoraxis)^2))

pdf("pred_parallel.pdf") 
plot(predmean(preds.multi), df.test$semimajoraxis)
dev.off()


round(c(rmse.default, rmse.P50, rmse.parallel),2)


###############################


#Prediction based on the best model () or the MPM (Median Probability Model)

get.best.model(result = result.default)
preds.best <- predict(get.best.model(result.default), df.test[, -1])
sqrt(mean((preds.best - df.test$semimajoraxis)^2))

get.mpm.model(result = result.default, y = df.train$semimajoraxis, x = df.train[, -1])
preds.mpm <- predict(get.mpm.model(result.default, y = df.train$semimajoraxis, x = df.train[, -1]), df.test[, -1])
sqrt(mean((preds.mpm - df.test$semimajoraxis)^2))



get.best.model(result = result.parallel)
preds.best_parallel <- predict(get.best.model(result.parallel), df.test[, -1])
sqrt(mean((preds.best_parallel - df.test$semimajoraxis)^2))



 
# Coefficients of the best model

coef(result.default)
coef(result.P50)
coef(result.parallel)


####################################################
#
# Diagnostic plots  (Section 3.8)
#
####################################################


pdf("diagn_default.pdf") 
diagn_plot(result.default, ylim = c(600,1500), FUN = max)
dev.off()
diagn_plot(result.default, ylim = c(600,1500), FUN = max)


pdf("diagn_long.pdf") 
diagn_plot(result.P50, ylim = c(600,1500), FUN = max)
dev.off()
diagn_plot(result.P50, ylim = c(600,1500), FUN = max)


pdf("diagn_par.pdf") 
diagn_plot(result.parallel, ylim = c(600,1500),FUN = max)
dev.off()

diagn_plot(result.parallel, ylim = c(600,1500),FUN = max)



#######################################################
#
# Example 2 (Section 4): 
#
# Zambia data set from the cAIC4 package
#
# Linear Mixed Model with Fractional Polynomials and other non-linear features
#
# Custom function to compute Marginal Likelihood with lme4, INLA and RTMB
#
#######################################################

rm(list = ls())

library(FBMS)

#install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
#options(repos=c( inlabruorg = "https://inlabru-org.r-universe.dev", INLA = "https://inla.r-inla-download.org/R/testing", CRAN = "https://cran.rstudio.com") )
#install.packages("fmesher") 

library(INLA)

library(tictoc)
library(lme4)
library(RTMB)

#install.packages("cAIC4") 
library(cAIC4)

data(Zambia, package = "cAIC4")
df <- as.data.frame(sapply(Zambia[1:5],scale))


transforms <- c("p0","p2","p3","p05","pm05","pm1","pm2","p0p0","p0p05","p0p1","p0p2","p0p3","p0p05","p0pm05","p0pm1","p0pm2")
probs <- gen.probs.gmjmcmc(transforms)
probs$gen <- c(1/3,1/3,0,1/3) # Modifications and interactions!

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


c(time.lme4$callback_msg, time.inla$callback_msg, time.rtmb$callback_msg)



#######################################################
#
# Serious analysis with lme4
#
#


# Analysis without non-linear features


result2a <- fbms(formula = z ~ 1+., data = df, N = 5000,
                 method = "mjmcmc.parallel", runs = 40, cores = parallel::detectCores()-1,
                 family = "custom", loglik.pi = mixed.model.loglik.lme4,
                 model_prior = list(r = 1/dim(df)[1]), 
                 extra_params = list(dr = droplevels(Zambia$dr)))

summary(result2a, labels = names(df)[-1])

plot(result2a)



# Analysis with fractional polynomials
probs <- gen.probs.gmjmcmc(transforms)
probs$gen <- c(1/3,1/3,0,1/3) # Modifications and interactions!

params <- gen.params.gmjmcmc(ncol(df) - 1)
params$feat$D <- 1   # Set depth of features to 1 (still allows for interactions)
params$feat$pop.max = 10

result2b <- fbms(formula = z ~ 1+., data = df, transforms = transforms,
                 probs = probs, params = params, P=25, N = 100,
                 method = "gmjmcmc.parallel", runs = 40, cores = parallel::detectCores()-1,
                 family = "custom", loglik.pi = mixed.model.loglik.lme4,
                 model_prior = list(r = 1/dim(df)[1]), 
                 extra_params = list(dr = droplevels(Zambia$dr)))

summary(result2b,tol = 0.05,labels=names(df)[-1])   



# Analysis with non-linear projections
transforms <- c("sigmoid")
probs <- gen.probs.gmjmcmc(transforms)
probs$gen <- c(0,0,0.5,0.5) 

params <- gen.params.gmjmcmc(ncol(df) - 1)
params$feat$pop.max = 10


result2c <- fbms(formula = z ~ 1+., data = df, transforms = transforms,
                 probs = probs, params = params, P=25, N = 100,
                 method = "gmjmcmc.parallel", runs = 40, cores = parallel::detectCores()-1,
                 family = "custom", loglik.pi = mixed.model.loglik.lme4,
                 model_prior = list(r = 1/dim(df)[1]), 
                 extra_params = list(dr = droplevels(Zambia$dr)))

summary(result2c,tol = 0.05,labels=names(df)[-1])   


# Comparison of results

summary(result2a, labels = names(df)[-1])
summary(result2b, labels = names(df)[-1])
summary(result2c, labels = names(df)[-1])


