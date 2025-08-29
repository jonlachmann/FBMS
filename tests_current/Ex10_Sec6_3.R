#######################################################
#
# Example 10 (Section 6.3): Epil data set from the INLA package
#
# Mixed Effect Poisson Model with Fractional Polynomials, using only fbms
#
# This is the valid version for the JSS Paper
#
#######################################################

#library(devtools)
#devtools::install_github("jonlachmann/FBMS@v1_arxiv", force=T, build_vignettes=F)


library(FBMS)
library(INLA)
library(tictoc)

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

# function to estimate log posterior
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

set.seed(03052024)
#specify indices for a random effect

result <- fbms(formula = y ~ 1+., data = df, transforms = transforms,
               method = "gmjmcmc", probs = probs, params = params, P=25, N = 100,
               family = "custom", loglik.pi = poisson.loglik.inla,
               model_prior = list(r = 1/dim(df)[1]), 
               extra_params = list(PID = data$Ind, INLA.num.threads = 1))

plot(result)
summary(result)


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