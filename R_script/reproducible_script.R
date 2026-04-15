###############################################################
# FBMS Reproducibility Script (JCS Submission)
# -------------------------------------------------------------
# This script reproduces examples in:
#
#   FBMS: Flexible Bayesian Model Selection and Model Averaging
#
# It installs the correct package versions and runs all the 
# examples used in the article.
#
# The script uses minimal, readable checks suitable for JCS:
#  - Mandatory packages are installed if missing
#  - Optional packages are installed if possible; otherwise skipped
#  - FBMS is always installed from a dedicated GitHub branch "JCS"
###############################################################



###############################################################
# 0. Install and load necessary packages
###############################################################

#  Ensure pak is available 
if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")

# Install/Update FBMS 

message("Installing FBMS and dependencies...")

packages_to_install <- c(
  "github::jonlachmann/FBMS@JCS", 
  "tictoc", 
  "kernlab",
  "lme4", 
  "cAIC4",
  "pec",
  "survival"
)

Sys.setenv(PAK_INPROCESS = "true") # Fixes the subprocess signal error
pak::pkg_install(packages_to_install, upgrade = TRUE)



library(FBMS)
library(tictoc)
library(parallel) # Base package, just needs loading

message("Environment ready for reproducibility.")

################################################################
################################################################
#
#  EXAMPLE 1: EXOPLANET DATA
#
#  Section 3 of the article
#
################################################################
################################################################

data(exoplanet)

train.indx <- 1:500
df.train = exoplanet[train.indx, ]
df.test  = exoplanet[-train.indx, ]

to3 <- function(x) x^3
transforms <- c("sigmoid", "sin_deg", "exp_dbl", "p0", "troot", "to3")


###############################################################
# Example 1.1 — Default single-thread GMJMCMC (Section 3)
###############################################################
set.seed(123)

result.default <- fbms(
  formula = semimajoraxis ~ 1 + .,
  data = df.train,
  method = "gmjmcmc",
  transforms = transforms
)


###############################################################
# Example 1.2 — Alternative priors 
###############################################################

set.seed(234)
result.BIC <- fbms(
  formula = semimajoraxis ~ 1 + .,
  data = df.train,
  method = "gmjmcmc",
  transforms = transforms,
  beta_prior = list(type = "Jeffreys-BIC", Var = "unknown")
)

set.seed(345)
result.EB <- fbms(
  formula = semimajoraxis ~ 1 + .,
  data = df.train,
  method = "gmjmcmc",
  transforms = transforms,
  beta_prior = list(type = "EB-global", a = 1)
)


###############################################################
# Example 1.3 — Longer single-thread run
###############################################################
set.seed(123)

result.P50 <- fbms(
  data = df.train,
  method = "gmjmcmc",
  transforms = transforms,
  P = 50, N = 1000, N.final = 5000
)


###############################################################
# Example 1.4 — Parallel GMJMCMC 
###############################################################
set.seed(123)

result.parallel <- fbms(
  data = df.train,
  method = "gmjmcmc.parallel",
  transforms = transforms,
  runs = 40,
  cores = parallel::detectCores() - 1,
  P = 25
)


###############################################################
# Example 1.5 — Summaries and plotting
###############################################################

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



plot(result.default)
plot(result.P50)
plot(result.parallel)

###############################################################
# Diagn_plot
# This part can be a bit time-consuming, so is by default not run
# If to run, change to diagnplot=TRUE
###############################################################
diagnplot=FALSE
if(diagnplot)
{
  set.seed(123)
  
  result.parallel2 <- fbms(
    data = df.train,
    method = "gmjmcmc.parallel",
    transforms = transforms,
    runs = 40,
    cores = parallel::detectCores() - 1,
    P = 100
  )
  
  diagn_plot(result.parallel2,type="min-mean-max",FUN=max)
  diagn_plot(result.parallel2,type="min-mean-max",mass=TRUE)
  diagn_plot(result.parallel2,type="total-mass")
}


###############################################################
# Example 1.6 — Prediction
###############################################################
preds <-  predict(result.default, df.test[,-1])
str(aggr(preds))

rmse.default <- sqrt(mean((predmean(preds) - df.test$semimajoraxis)^2))
plot(predmean(preds), df.test$semimajoraxis)


###############################

preds.P50 <- predict(result.P50, df.test[,-1])  
rmse.P50 <-  sqrt(mean((predmean(preds.P50) - df.test$semimajoraxis)^2))
plot(predmean(preds.P50), df.test$semimajoraxis)


###############################


preds.multi <- predict(result.parallel , df.test[,-1], link = function(x) x)
rmse.parallel <- sqrt(mean((predmean(preds.multi) - df.test$semimajoraxis)^2))
plot(predmean(preds.multi), df.test$semimajoraxis)


round(c(rmse.default, rmse.P50, rmse.parallel),2)


###############################################################
# Example 1.7 — Best model & MPM
###############################################################

best.default <- get.best.model(result.default)
mpm.default  <- get.mpm.model(result.default,
                              y = df.train$semimajoraxis,
                              x = df.train[,-1])

sqrt(mean((predict(best.default, df.test[,-1]) -
             df.test$semimajoraxis)^2))

sqrt(mean((predict(mpm.default, df.test[,-1]) -
             df.test$semimajoraxis)^2))





################################################################
################################################################
#
#  EXAMPLE 2: Logistic Regression
#
#  Section 4.1 of the article
#
################################################################
################################################################

rm(list = ls())
library(FBMS)
library(kernlab)


###############################################################
# 2.0 Load spam data (requires kernlab)
###############################################################

data("spam")
df <- spam[,c(58,1:57)]

#number of observations and covariates

n <- dim(df)[1] 
p <- dim(df)[2] - 1   

colnames(df) <-  c("y", paste0("x",1:p))
df$y = as.numeric(df$y == "spam")

to3 <- function(x) x^3
transforms <- c("sigmoid","sin_deg","exp_dbl","p0","troot","to3")
probs <- gen.probs.gmjmcmc(transforms)
probs$gen <- c(1,1,1,1) 

params <- gen.params.gmjmcmc(p)
params$feat$check.col <- F

####################################################
# 2.1 Single thread analysis
####################################################


set.seed(6001)
# Perform analysis with logistic.loglik
result <- fbms(formula = y~1+.,data = df, method = "gmjmcmc", 
               family = "binomial", beta_prior = list(type = "Jeffreys-BIC"),
               transforms = transforms, probs = probs, params = params)

summary(result)


#######################
#
# Prediction accuracy
# IMPORTANT: specify correct link function for predict
#

# Model averaging
pred <- predict(result, x =  df[,-1], link = function(x)(1/(1+exp(-x))))  
mean(round(pred$aggr$mean)==df$y)

# Best model
bm <- get.best.model(result = result)
preds <-  predict(object = bm, df[,-1],link = function(x)(1/(1+exp(-x))))
mean(round(preds)==df$y)

# Median Probability Model
mpm <- get.mpm.model(result = result,family = "binomial",y = df$y,x=df[,-1])
preds <-  predict(mpm, df[,-1],link = function(x)(1/(1+exp(-x))))
mean(round(preds)==df$y)



plot(pred$aggr$mean)
points(pred$aggr$quantiles[1,], col = 2)
points(pred$aggr$quantiles[3,], col = 3)


head(cbind(pred$aggr$mean, pred$aggr$quantiles[1,],pred$aggr$quantiles[3,]))


####################################################
# 2.2 Multiple thread analysis
####################################################

set.seed(6002)

result_parallel <- fbms(formula = y~1+.,data = df, method = "gmjmcmc.parallel", 
                        family = "binomial", beta_prior = list(type = "Jeffreys-BIC"),
                        runs = 40, cores = parallel::detectCores() - 1, transforms = transforms, 
                        probs = probs, params = params, P=25)

summary(result_parallel)

#######################
#
# Prediction accuracy
# IMPORTANT: specify correct link function for predict
#

# Model averaging
pred_parallel = predict(result_parallel, x =  df[,-1], link = function(x)(1/(1+exp(-x))))  
mean(round(pred_parallel$aggr$mean)==df$y)

# Best Model
#bm_parallel <- get.best.model(result_parallel)
#pred_bm_parallel <-  predict(bm_parallel, df[,-1],link = function(x)(1/(1+exp(-x))))
#mean(round(pred_bm_parallel)==df$y)

# Median Probability Model
mpm_parallel <-  predict(get.mpm.model(result = result_parallel,family = "binomial",y = df$y,x=df[,-1]), df[,-1],link = function(x)(1/(1+exp(-x))))
mean(round(mpm_parallel)==df$y)



####################################################
# multiple thread analysis - random strategy (fully non-linear)
####################################################

set.seed(6003)

params$feat$alpha = "random"

result_parallel_random <- fbms(formula = y~1+.,data = df, method = "gmjmcmc.parallel", 
                               family = "binomial", beta_prior = list(type = "Jeffreys-BIC"),
                               runs = 40, cores = parallel::detectCores() - 1, transforms = transforms, 
                               probs = probs, params = params, P=25)

summary(result_parallel_random)

#######################
#
# Prediction accuracy
# IMPORTANT: specify correct link function for predict
#

# Model averaging
pred_parallel_random = predict(result_parallel_random, x =  df[,-1], link = function(x)(1/(1+exp(-x))))  
mean(round(pred_parallel_random$aggr$mean)==df$y)

# Median Probability Model
mpm_parallel_random <-  predict(get.mpm.model(result = result_parallel_random,family = "binomial",y = df$y,x=df[,-1]), df[,-1],link = function(x)(1/(1+exp(-x))))
mean(round(mpm_parallel_random)==df$y)






################################################################
################################################################
#
#  EXAMPLE 3: MIXED MODELS WITH FRACTIONAL POLYNOMIALS
#
#  Section 4.2 of the article
#
################################################################
################################################################

rm(list = ls())
library(FBMS)


###############################################################
# 3.0 Load Zambia data (requires cAIC4)
###############################################################
if (!requireNamespace("lme4", quietly = TRUE)) {
  stop("Optional package 'lme4' is required for Example 2. Please install it.")
}
if (!requireNamespace("tictoc", quietly = TRUE)) {
  stop("Optional package 'tictoc' is required for Example 2. Please install it.")
}
if (!requireNamespace("cAIC4", quietly = TRUE)) {
  stop("Optional package 'cAIC4' is required for Example 2. Please install it.")
}

library(tictoc)
library(lme4)


data(Zambia, package = "cAIC4")
df <- as.data.frame(sapply(Zambia[1:5],scale))


transforms <- c("p0","p2","p3","p05","pm05","pm1","pm2",
                "p0p0","p0p05","p0p1","p0p2","p0p3",
                "p0p05","p0pm05","p0pm1","p0pm2")


probs <- gen.probs.gmjmcmc(transforms)
probs$gen <- c(1/3,1/3,0,1/3) # Modifications and interactions!

params <- gen.params.gmjmcmc(ncol(df) - 1)
params$feat$D <- 1   # Set depth of features to 1 (still allows for interactions)
params$feat$pop.max <- 10



###############################################################
# 3.1 Define custom log-likelihoods for lme4
###############################################################

# lme4 version
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
  
  mloglik <- as.numeric(logLik(mm))  -  0.5*log(length(y)) * (dim(data)[2] - 2) #Laplace approximation
  
  # logarithm of model prior
  if (length(mlpost_params$r) == 0)  mlpost_params$r <- 1/dim(x)[1]  # default value or parameter r
  lp <- log_prior(mlpost_params, complex)
  
  
  return(list(crit = mloglik + lp, coefs = fixef(mm)))
}

###############################################################
# 3.2 Small demonstration run for runtime comparisons
###############################################################


set.seed(03052024)

tic()
result1a <- fbms(formula = z ~ 1+., data = df, transforms = transforms,
                 method = "gmjmcmc",probs = probs, params = params, P=3, N = 30,
                 family = "custom", loglik.pi = mixed.model.loglik.lme4,
                 model_prior = list(r = 1/dim(df)[1]), 
                 extra_params = list(dr = droplevels(Zambia$dr)))
time.lme4 <- toc()


cat(c(time.lme4$callback_msg))

###############################################################
# 3.3 Serious analysis with lme4 (Section 4). Runs within time
# constraints of JCS on Apple M1 Max 32 GB, but can be slower
# on older machines. Please, set run.long.mixed = FALSE, if the
# example exceeds reasonable time.
###############################################################

# Specify if to run long chains under mixed effect models.
# Default is false as these chains an run longer than 20 minutes 
# depending on the machines used. 
run.long.mixed <- TRUE

if(run.long.mixed)
{
  probs <- gen.probs.gmjmcmc(transforms)
  params <- gen.params.gmjmcmc(ncol(df) - 1)
  params$feat$D <- 1
  params$feat$pop.max <- 10
  
  
  # No nonlinear features
  result2a <- fbms(
    formula = z ~ 1+., data = df,
    N = 5000,
    method = "mjmcmc.parallel",
    runs = 40,
    cores = parallel::detectCores() - 1,
    family = "custom",
    loglik.pi = mixed.model.loglik.lme4,
    model_prior = list(r = 1/nrow(df)),
    extra_params = list(dr = droplevels(Zambia$dr))
  )
  
  summary(result2a, labels = names(df)[-1])
  plot(result2a)
  
  
  # Fractional polynomials
  result2b <- fbms(
    formula = z ~ 1+., data = df,
    transforms = transforms, probs = probs, params = params,
    P = 25, N = 100,
    method = "gmjmcmc.parallel",
    runs = 40,
    cores = parallel::detectCores() - 1,
    family = "custom",
    loglik.pi = mixed.model.loglik.lme4,
    model_prior = list(r = 1/nrow(df)),
    extra_params = list(dr = droplevels(Zambia$dr))
  )
  
  summary(result2b, tol = 0.05, labels = names(df)[-1])
  
  
  # Non-linear projections
  transforms.sigmoid <- c("sigmoid")
  probs.sigmoid <- gen.probs.gmjmcmc(transforms.sigmoid)
  probs.sigmoid$gen <- c(0, 0, 0.5, 0.5)
  
  params.sigmoid <- gen.params.gmjmcmc(ncol(df) - 1)
  params.sigmoid$feat$pop.max <- 10
  
  result2c <- fbms(
    formula = z ~ 1+., data = df,
    transforms = transforms.sigmoid,
    probs = probs.sigmoid,
    params = params.sigmoid,
    P = 25, N = 100,
    method = "gmjmcmc.parallel",
    runs = 40,
    cores = parallel::detectCores() - 1,
    family = "custom",
    loglik.pi = mixed.model.loglik.lme4,
    model_prior = list(r = 1/nrow(df)),
    extra_params = list(dr = droplevels(Zambia$dr))
  )
  
  summary(result2c, tol = 0.05, labels = names(df)[-1])
  
  
  # Comparison
  summary(result2a, labels = names(df)[-1])
  summary(result2b, labels = names(df)[-1])
  summary(result2c, labels = names(df)[-1])
}


################################################################
################################################################
#
#  EXAMPLE 5: COX REGRESSION
#
#  Section 4.4 of the article
#
################################################################
################################################################


rm(list = ls())

library(FBMS)
library(pec) #for the computation of cindex
library(survival)



###############################################################
# 5.0 Load GSGB breast cancer dataset and data preparation
###############################################################

download.file('https://www.uniklinik-freiburg.de/fileadmin/mediapool/08_institute/biometrie-statistik/Dateien/Studium_und_Lehre/Lehrbuecher/Multivariable_Model-building/gbsg_br_ca.zip',
              'gbsg_br_ca.zip')
df1 <- read.csv(unz('gbsg_br_ca.zip',
                    'gbsg_br_ca/gbsg_br_ca.csv'),
                header = TRUE)

# Prepare data
df <- df1[, c(13, 14, 2:4, 6:8, 10:12)]
names(df) = c("time","cens",names(df)[3:ncol(df)])

# Split into training and test set
set.seed(123)
train <- c(sample((1:nrow(df))[df$cens == 1], sum(df$cens)*2/3), # split separately events
           sample((1:nrow(df))[df$cens == 0], sum(!df$cens)*2/3)) # and censored observations

df.train <- df[train,]
df.test <- df[-train,]

# time will be used as an extra parameter in the custom function
time <- df.train$time


params <- gen.params.gmjmcmc(ncol(df.train) - 2)
transforms <- c("p0","p2","p3","p05","pm05","pm1","pm2")
probs <- gen.probs.gmjmcmc(transforms)


###############################################################
# 5.1 Define custom log-likelihoods for Cox regression
###############################################################

surv.pseudo.loglik = function(y, x, model, complex, mlpost_params){
  
  data <- data.frame(time = mlpost_params$time, cens = y, as.matrix(x[,model]))[,-3]  # Removing intercept
  if(dim(data)[2]==2)
  {  
    formula1 <- as.formula(paste0("Surv(time,cens)","~ 1"))
    out <- coxph(formula1, data = data)
    out$loglik <- c(out$loglik,out$loglik)
    out$coefficients <- NULL
  }  else {
    formula1 <- as.formula(paste0("Surv(time,cens)","~ 1 + ."))
    
    out <- coxph(formula1, data = data)
  }  
  # logarithm of marginal likelihood
  mloglik <- (out$loglik[2] - out$loglik[1]) -  log(length(y)) * (dim(data)[2] - 2)/2   
  
  # logarithm of model prior
  if (length(mlpost_params$r) == 0)  mlpost_params$r <- 1/dim(x)[1]  # default value or parameter r
  lp <- log_prior(mlpost_params, complex)
  
  # Compute criterion and consider special cases like multicollinearity
  
  crit <- mloglik + lp
  if(sum(is.na(out$coefficients))>0)   #Get rid of models with collinearities (with more than two features)
    crit <- -.Machine$double.xmax
  
  return(list(crit = crit, coefs =  c(0,out$coefficients)))
  
}


###############################################################
# 5.2 Analysis with 4 different modeling strategies
###############################################################


# 1) Single chain analysis (just to illustrate how it works)
set.seed(121)
result1 <- fbms(formula = cens ~ 1 + .,data = df.train[,-1], P = 5,
                transforms = transforms, method = "gmjmcmc",
                family = "custom", loglik.pi = surv.pseudo.loglik, 
                model_prior = list(r = 0.5),
                extra_params = list(time = time))


summary(result1, data = df, tol = 0.01, effects = c(0.025,0.5,0.975))

# 2) Parallel version only linear terms with mjmcmc


set.seed(122)

result2 <- fbms(formula = cens ~ 1 + .,data = df.train[,-1], params = params,  
                method = "mjmcmc.parallel", runs = 40,  cores = parallel::detectCores() - 1,
                family = "custom", loglik.pi = surv.pseudo.loglik, 
                model_prior = list(r = 0.5), extra_params = list(time = time))

summary(result2,tol = 0.01,labels = names(df.train)[-(1:2)],effects = c(0.025,0.5,0.975))



# 3) Parallel version only fractional polynomials

set.seed(123)
probs$gen <- c(0,1,0,1)
params$feat$D = 1

result3 <- fbms(formula = cens ~ 1 + .,data = df.train[,-1], params = params, probs = probs, P = 10,
                transforms = transforms, method = "gmjmcmc.parallel", runs = 40,  cores = parallel::detectCores() - 1,
                family = "custom", loglik.pi = surv.pseudo.loglik, 
                model_prior = list(r = 0.5), extra_params = list(time = time))


summary(result3,tol = 0.01, effects = c(0.025,0.5,0.975))



# 4) Parallel version using all types of non-linear features
set.seed(124)
probs$gen <- c(1,1,1,1)
params$feat$D = 5

result4 <- fbms(formula = cens ~ 1 + .,data = df.train[,-1], params = params, probs = probs,P = 10, 
                transforms = transforms, method = "gmjmcmc.parallel", runs = 40,  cores = parallel::detectCores() - 1, N = 300,
                family = "custom", loglik.pi = surv.pseudo.loglik, 
                model_prior = list(r = 0.5), extra_params = list(time = time))


summary(result4,tol = 0.01,effects = c(0.025,0.5,0.975))



####################################################
#  5.3 Prediction and C index using model averaging
####################################################


linpreds1.train <- predict(result1,df.train[,-(1:2)], link = function(x) x)
linpreds1 <- predict(result1,df.test[,-(1:2)], link = function(x) x)

linpreds2.train <- predict(result2,df.train[,-(1:2)], link = function(x) x)
linpreds2 <- predict(result2,df.test[,-(1:2)], link = function(x) x)

linpreds3.train <- predict(result3,df.train[,-(1:2)], link = function(x) x)
linpreds3 <- predict(result3,df.test[,-(1:2)], link = function(x) x)

linpreds4.train <- predict(result4,df.train[,-(1:2)], link = function(x) x)
linpreds4 <- predict(result4,df.test[,-(1:2)], link = function(x) x)



df.train$average.lin.pred1 <- linpreds1.train$aggr$mean
df.train$average.lin.pred2 <- linpreds2.train$aggr$mean
df.train$average.lin.pred3 <- linpreds3.train$aggr$mean
df.train$average.lin.pred4 <- linpreds4.train$aggr$mean

df.test$average.lin.pred1 <- linpreds1$aggr$mean
df.test$average.lin.pred2 <- linpreds2$aggr$mean
df.test$average.lin.pred3 <- linpreds3$aggr$mean
df.test$average.lin.pred4 <- linpreds4$aggr$mean



# Compute cindex using package pec

mod1 <- coxph(Surv(time, cens) ~ average.lin.pred1, data = as.data.frame(df.train), x = TRUE)
cindex1 <- cindex(mod1, mod1$formula, data = as.data.frame(df.test), cens.model = 'cox')$AppCindex

mod2 <- coxph(Surv(time, cens) ~ average.lin.pred2, data = as.data.frame(df.train), x = TRUE)
cindex2 <- cindex(mod2, mod2$formula, data = as.data.frame(df.test), cens.model = 'cox')$AppCindex

mod3 <- coxph(Surv(time, cens) ~ average.lin.pred3, data = as.data.frame(df.train), x = TRUE)
cindex3 <- cindex(mod3, mod3$formula, data = as.data.frame(df.test), cens.model = 'cox')$AppCindex

mod4 <- coxph(Surv(time, cens) ~ average.lin.pred4, data = as.data.frame(df.train),  x = TRUE)
cindex4 <- cindex(mod4, mod4$formula, data = as.data.frame(df.test), cens.model = 'cox')$AppCindex


#Full model without nonlinearities (for the sake of comparison)
mod5 <- coxph(Surv(time, cens) ~ 1+., data = as.data.frame(df.train[,1:11]),x = T)
cindex5 <- cindex(mod5, mod5$formula, data = as.data.frame(df.test), cens.model = 'cox')$AppCindex

#Model without predictors (for the sake of comparison)
mod6 <- coxph(Surv(time, cens) ~ 1, data = as.data.frame(df.train[,1:11]),x = T)
cindex6 <- cindex(mod6, mod6$formula, data = as.data.frame(df.test), cens.model = 'cox')$AppCindex

all.cindices = round(unlist(c(cindex1, cindex2, cindex3, cindex4, cindex5, cindex6)),3)
names(all.cindices) = c("Model 1", "Model 2", "Model 3", "Model 4", "Full Linear Model", "Null Model")



####################################################
#  5.3 Prediction and C index using best model 
####################################################


linpreds.train.best <- predict(get.best.model(result1),df.train[,-(1:2)], link = function(x) x)
linpreds.best <- predict(get.best.model(result1),df.test[,-(1:2)], link = function(x) x)


linpreds2.train.best <- predict(get.best.model(result2),df.train[,-(1:2)], link = function(x) x)
linpreds2.best <- predict(get.best.model(result2),df.test[,-(1:2)], link = function(x) x)


linpreds3.train.best <- predict(get.best.model(result3),df.train[,-(1:2)], link = function(x) x)
linpreds3.best <- predict(get.best.model(result3),df.test[,-(1:2)], link = function(x) x)


linpreds4.train.best <- predict(get.best.model(result4),df.train[,-(1:2)], link = function(x) x)
linpreds4.best <- predict(get.best.model(result4),df.test[,-(1:2)], link = function(x) x)


df.train$best.lin.pred1 <- linpreds.train.best
df.train$best.lin.pred2 <- linpreds2.train.best
df.train$best.lin.pred3 <- linpreds3.train.best
df.train$best.lin.pred4 <- linpreds4.train.best

df.test$best.lin.pred1 <- linpreds.best
df.test$best.lin.pred2 <- linpreds2.best
df.test$best.lin.pred3 <- linpreds3.best
df.test$best.lin.pred4 <- linpreds4.best

mod1 <- coxph(Surv(time, cens) ~ best.lin.pred1, data = as.data.frame(df.train), x = TRUE)
cindex1 <- cindex(mod1, mod1$formula, data = as.data.frame(df.test), cens.model = 'cox')$AppCindex

mod2 <- coxph(Surv(time, cens) ~ best.lin.pred2, data = as.data.frame(df.train), x = TRUE)
cindex2 <- cindex(mod2, mod2$formula, data = as.data.frame(df.test), cens.model = 'cox')$AppCindex

mod3 <- coxph(Surv(time, cens) ~ best.lin.pred3, data = as.data.frame(df.train), x = TRUE)
cindex3 <- cindex(mod3, mod3$formula, data = as.data.frame(df.test), cens.model = 'cox')$AppCindex

mod4 <- coxph(Surv(time, cens) ~ best.lin.pred4, data = as.data.frame(df.train),  x = TRUE)
cindex4 <- cindex(mod4, mod4$formula, data = as.data.frame(df.test), cens.model = 'cox')$AppCindex



all.cindices = rbind(all.cindices, round(unlist(c(cindex1, cindex2, cindex3, cindex4, cindex5, cindex6)),3))


####################################################
#  5.4 Prediction and C index using mpm model 
####################################################


linpreds.train.mpm <- predict(get.mpm.model(result1, y = df.train$cens, x = df.train[, -c(1,2)],family = "custom",
                                            loglik.pi = surv.pseudo.loglik,params = list(r = 0.5, time = time)),
                              df.train[,-(1:2)], link = function(x) x)

linpreds.mpm <- predict(get.mpm.model(result1, y = df.train$cens, x = df.train[, -c(1,2)],family = "custom",
                                      loglik.pi = surv.pseudo.loglik,params = list(r = 0.5, time = time)),df.test[,-(1:2)], link = function(x) x)


linpreds2.train.mpm <- predict(get.mpm.model(result2, y = df.train$cens, x = df.train[, -c(1,2)],family = "custom",
                                             loglik.pi = surv.pseudo.loglik,params = list(r = 0.5, time = time)),
                               df.train[,-(1:2)], link = function(x) x)
linpreds2.mpm <- predict(get.mpm.model(result2, y = df.train$cens, x = df.train[, -c(1,2)],family = "custom",
                                       loglik.pi = surv.pseudo.loglik,params = list(r = 0.5, time = time)),df.test[,-(1:2)], link = function(x) x)

linpreds3.train.mpm <- predict(get.mpm.model(result3, y = df.train$cens, x = df.train[, -c(1,2)],family = "custom",
                                             loglik.pi = surv.pseudo.loglik,params = list(r = 0.5, time = time)),
                               df.train[,-(1:2)], link = function(x) x)
linpreds3.mpm <- predict(get.mpm.model(result3, y = df.train$cens, x = df.train[, -c(1,2)],family = "custom",
                                       loglik.pi = surv.pseudo.loglik,params = list(r = 0.5, time = time)),df.test[,-(1:2)], link = function(x) x)


linpreds4.train.mpm <- predict(get.mpm.model(result4, y = df.train$cens, x = df.train[, -c(1,2)],family = "custom",
                                             loglik.pi = surv.pseudo.loglik,params = list(r = 0.5, time = time)),
                               df.train[,-(1:2)], link = function(x) x)
linpreds4.mpm <- predict(get.mpm.model(result4, y = df.train$cens, x = df.train[, -c(1,2)],family = "custom",
                                       loglik.pi = surv.pseudo.loglik,params = list(r = 0.5, time = time)),df.test[,-(1:2)], link = function(x) x)


df.train$mpm.lin.pred1 <- linpreds.train.mpm
df.train$mpm.lin.pred2 <- linpreds2.train.mpm
df.train$mpm.lin.pred3 <- linpreds3.train.mpm
df.train$mpm.lin.pred4 <- linpreds4.train.mpm

df.test$mpm.lin.pred1 <- linpreds.mpm
df.test$mpm.lin.pred2 <- linpreds2.mpm
df.test$mpm.lin.pred3 <- linpreds3.mpm
df.test$mpm.lin.pred4 <- linpreds4.mpm

mod1 <- coxph(Surv(time, cens) ~ mpm.lin.pred1, data = as.data.frame(df.train), x = TRUE)
cindex1 <- cindex(mod1, mod1$formula, data = as.data.frame(df.test), cens.model = 'cox')$AppCindex

mod2 <- coxph(Surv(time, cens) ~ mpm.lin.pred2, data = as.data.frame(df.train), x = TRUE)
cindex2 <- cindex(mod2, mod2$formula, data = as.data.frame(df.test), cens.model = 'cox')$AppCindex

mod3 <- coxph(Surv(time, cens) ~ mpm.lin.pred3, data = as.data.frame(df.train), x = TRUE)
cindex3 <- cindex(mod3, mod3$formula, data = as.data.frame(df.test), cens.model = 'cox')$AppCindex

mod4 <- coxph(Surv(time, cens) ~ mpm.lin.pred4, data = as.data.frame(df.train),  x = TRUE)
cindex4 <- cindex(mod4, mod4$formula, data = as.data.frame(df.test), cens.model = 'cox')$AppCindex


all.cindices = rbind(all.cindices, round(unlist(c(cindex1, cindex2, cindex3, cindex4, cindex5, cindex6)),3))
rownames(all.cindices) = c("Model Averaging", "Best Model", "MPM")

all.cindices


