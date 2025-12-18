###############################################################
# FBMS Reproducibility Script (JSS Submission)
# -------------------------------------------------------------
# This script reproduces examples in:
#
#   FBMS: Flexible Bayesian Model Selection and Model Averaging
#
# It installs the correct package versions and runs the two
# main examples used in the article.
#
# The script uses minimal, readable checks suitable for JSS:
#  - Mandatory packages are installed if missing
#  - Optional packages are installed if possible; otherwise skipped
#  - FBMS is always installed from a dedicated GitHub branch "jss_v2"
###############################################################



###############################################################
# 0. Helper: Install a mandatory package (stop if fails)
###############################################################

install_mandatory <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing mandatory package: ", pkg)
    tryCatch(
      install.packages(pkg),
      error = function(e) {
        stop("Failed to install mandatory package ", pkg, call. = FALSE)
      }
    )
  }
}

###############################################################
# 1. Install mandatory packages
###############################################################

mandatory_pkgs <- c("devtools", "parallel", "tictoc", "lme4","cAIC4")

for (p in mandatory_pkgs) install_mandatory(p)

library(devtools)


###############################################################
# 2. Install FBMS (always from GitHub to enforce correct version)
###############################################################

message("Installing FBMS from GitHub (branch softwareX)...")
install_github("jonlachmann/FBMS@jsoftwareX",
               force = TRUE, build_vignettes = FALSE)

library(FBMS)
library(tictoc)

################################################################
################################################################
#
#  EXAMPLE 1: EXOPLANET DATA
#
#  Section 3 of the article
#
################################################################
################################################################

library(FBMS)
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
#  EXAMPLE 2: MIXED MODELS WITH FRACTIONAL POLYNOMIALS
#
#  Section 4 of the article
#
################################################################
################################################################

rm(list = ls())
library(FBMS)


###############################################################
# 2.0 Load Zambia data (requires cAIC4)
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
# 2.1 Define custom log-likelihoods for lme4
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
  
  mloglik <- as.numeric(logLik(mm))  -  0.5*log(length(y)) * (dim(data)[2] - 2) #Laplace approximation for beta prior
  
  # logarithm of model prior
  if (length(mlpost_params$r) == 0)  mlpost_params$r <- 1/dim(x)[1]  # default value or parameter r
  lp <- log_prior(mlpost_params, complex)
  
  
  return(list(crit = mloglik + lp, coefs = fixef(mm)))
}

###############################################################
# 2.2 Small demonstration run for runtime comparisons
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
# 2.3 Serious analysis with lme4 (Section 4). Runs within time
# constraints of JSS on Apple M1 Max 32 GB, but can be slower
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

