
####################################################################
# Install INLA and RTMB packages (consult with IT if installation fails, 
# which may occasionally happen for INLA as it is not on CRAN). JCS branch
###################################################################



if (!requireNamespace("RTMB", quietly = TRUE)) {
  message("Trying to install optional package RTMB...")
  try(install.packages("RTMB"), silent = TRUE)
}

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
#  EXAMPLE 2: MIXED MODELS WITH FRACTIONAL POLYNOMIALS
#
#  Section 4 of the article
#
################################################################
################################################################


library(FBMS)



###############################################################
# 2.0 Load Zambia data (requires cAIC4)
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
# 2.1 Define custom log-likelihoods for INLA, RTMB
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
# 2.2 Small demonstration run for runtime comparisons
###############################################################

set.seed(3052024)

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

