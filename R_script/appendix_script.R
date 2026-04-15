
###############################################################
# 0. Install and load necessary packages
###############################################################


if (!requireNamespace("pak", quietly = TRUE)) install.packages("pak")

# Install/Update FBMS and other standard packages needed

message("Installing FBMS and dependencies...")

packages_to_install <- c(
  "github::jonlachmann/FBMS@JCS", 
  "tictoc", 
  "lme4", 
  "cAIC4"
)

pak::pkg_install(packages_to_install, upgrade = TRUE)




# Install RTMB (Standard CRAN package)
if (!requireNamespace("RTMB", quietly = TRUE)) {
  message("Attempting to install optional package: RTMB")
  try(pak::pkg_install("RTMB"), silent = TRUE)
}


####################################################################
# Install INLA packages (consult with IT if installation fails, 
# which may occasionally happen for INLA as it is not on CRAN). 
###################################################################


if (!requireNamespace("INLA", quietly = TRUE)) {
  message("Trying to install optional package INLA...")
  
  # Try to load the installer (only if previously installed)
  if (!requireNamespace("INLA", quietly = TRUE)) {
    tryCatch(
      {
        install.packages("INLA",repos=c(getOption("repos"),INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
      },
      error = function(e) {
        message("INLA could not be installed; continuing without it.")
      }
    )
  }
}



################################################################
################################################################
#
#  EXAMPLE 3A: LINEAR MIXED MODELS WITH FRACTIONAL POLYNOMIALS
#
#  Analysis with INLA and RTMB
#
#  Section 4.2.2 of the article
#
################################################################
################################################################


library(FBMS)



###############################################################
# 3A.0 Load Zambia data (requires cAIC4)
###############################################################
if (!requireNamespace("cAIC4", quietly = TRUE)) {
  stop("Optional package 'cAIC4' is required for Example 2. Please install it.")
}

data(Zambia, package = "cAIC4")
df <- as.data.frame(sapply(Zambia[1:5],scale))


transforms <- c("p0","p2","p3","p05","pm05","pm1","pm2",
                "p0p0","p0p05","p0p1","p0p2","p0p3",
                "p0p05","p0pm05","p0pm1","p0pm2")


probs <- gen.probs.gmjmcmc(transforms)
probs$gen <- c(1/3,1/3,0,1/3) # Modifications and interactions!

params <- gen.params.gmjmcmc(ncol(df) - 1)
params$feat$D <- 1   # Set depth of features to 1 (still allows for interactions)
params$feat$pop.max = 10



###############################################################
# 3A.1 Define custom log-likelihoods for INLA, RTMB
###############################################################

# ---------------------------------------------------------------
# INLA version (only used if INLA is properly installed)
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

# ---------------------------------------------------------------
# RTMB version (only used if RTMB is properly installed)
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


###############################################################
# 3A.2 Small demonstration run for runtime comparisons
###############################################################

set.seed(3052026)

library(tictoc)


time.inla <- -1
if (requireNamespace("INLA", quietly = TRUE)) {
  library(INLA)
  library(cAIC4)
  
  data(Zambia, package = "cAIC4")
  df <- as.data.frame(sapply(Zambia[1:5],scale))
  
  tic()
  result1b <- fbms(
    formula = z ~ 1+., data = df,
    transforms = transforms,
    method = "gmjmcmc", P = 3, N = 30,
    family = "custom",
    loglik.pi = mixed.model.loglik.inla,
    model_prior = list(r = 1/nrow(df)),
    extra_params = list(dr = droplevels(Zambia$dr),
                        INLA.num.threads = 4)
  )
  time.inla <- toc()
}

time.rtmb <- -1
if (requireNamespace("RTMB", quietly = TRUE)) {
  library(RTMB)
  
  data(Zambia, package = "cAIC4")
  df <- as.data.frame(sapply(Zambia[1:5],scale))
  
  
  tic()
  result1c <- fbms(
    formula = z ~ 1+., data = df,
    transforms = transforms,
    method = "gmjmcmc", P = 3, N = 30,
    family = "custom",
    loglik.pi = mixed.model.loglik.rtmb,
    model_prior = list(r = 1/nrow(df)),
    extra_params = list(
      dr = droplevels(Zambia$dr),
      nr_dr = sum(table(Zambia$dr) > 0)
    )
  )
  time.rtmb <- toc()
}

cat(c(time.inla$callback_msg, time.rtmb$callback_msg))



################################################################
################################################################
#
#  EXAMPLE 4: MIXED EFFECT POISSON MODEL
#
#  Section 4.3 of the article
#
################################################################
################################################################

rm(list = ls())

library(FBMS)
library(INLA)
library(tictoc)

###############################################################
# 4.0 Load Epil data (requires INLA)
###############################################################


data <- INLA::Epil
data <- data[,-c(5,6)]

df <- data[1:5]
df$V2 <- rep(c(0,1,0,0),59)
df$V3 <- rep(c(0,0,1,0),59)
df$V4 <- rep(c(0,0,0,1),59)

transforms <- c("p0","p2","p3","p05","pm05","pm1","pm2","p0p0","p0p05","p0p1","p0p2","p0p3","p0p05","p0pm05","p0pm1","p0pm2")
probs <- gen.probs.gmjmcmc(transforms)
probs$gen <- c(1,1,0,1) # Only modifications!

params <- gen.params.gmjmcmc(ncol(df) - 1)
params$feat$D <- 2   # Set depth of features to 2 (allow for interactions)
params$feat$keep.min <- 0.2
params$greedy$steps <- 2
params$greedy$tries <- 1
params$sa$t.min <- 0.1
params$sa$dt <- 10


###############################################################
# 4.1 Define custom log-likelihoods for INLA 
###############################################################

poisson.loglik.inla <- function (y, x, model, complex, mlpost_params) 
{
  if(sum(model)>1)
  {
    data1 <- data.frame(y, as.matrix(x[,model]), mlpost_params$PID)
    formula1 <- as.formula(paste0(names(data1)[1],"~",paste0(names(data1)[3:(dim(data1)[2]-1)],collapse = "+"),"+ f(mlpost_params.PID,model = \"iid\")"))
  } else
  {
    data1 <- data.frame(y, mlpost_params$PID)
    formula1 <- as.formula(paste0(names(data1)[1],"~","1 + f(mlpost_params.PID,model = \"iid\")"))
  }
  
  #to make sure inla is not stuck
  inla.setOption(inla.timeout=30)
  inla.setOption(num.threads=mlpost_params$INLA.num.threads) 
  
  mod<-NULL
  
  #error handling for unstable libraries that might crash
  tryCatch({
    mod <- inla(family = "poisson",silent = 1L,safe = F,data = data1,formula = formula1)
  }, error = function(e) {
    # Handle the error by setting result to NULL
    mod <- NULL
    # Print a message or log the error if needed
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



###############################################################
# 4.2 Fast gmjmcmc run with P = 3
###############################################################


set.seed(03052024)
#specify indices for a random effect

result <- fbms(formula = y ~ 1+., data = df, transforms = transforms,
               method = "gmjmcmc", probs = probs, params = params, P=3, N = 100,
               family = "custom", loglik.pi = poisson.loglik.inla,
               model_prior = list(r = 1/dim(df)[1]), 
               extra_params = list(PID = data$Ind, INLA.num.threads = 1))

plot(result)
summary(result)


###############################################################
# 4.3 Full analysis with gmjmcmc.parallel
###############################################################


set.seed(23052024)

tic()
# Number of threads used by INLA set to 1 to avoid conflicts between two layers of parallelization
result2 <- fbms(formula = y ~ 1+., data = df, transforms = transforms,
                probs = probs, params = params, P=25, N = 100,
                method = "gmjmcmc.parallel", runs = 40, cores = 40,
                family = "custom", loglik.pi = poisson.loglik.inla,
                model_prior = list(r = 1/dim(df)[1]), 
                extra_params = list(PID = data$Ind, INLA.num.threads = 1))
time.inla <- toc()

plot(result2)
summary(result2, labels = names(df)[-1], tol = 0.01)

