#######################################################
#
# Example 11 (Section 6.4):
#
# Subsampling, using only fbms
#
# Heart Disease Health Indicators Dataset”
#
# This is the valid version for the JSS Paper
#
#######################################################

#library(devtools)
#devtools::install_github("jonlachmann/FBMS@v1_arxiv", force=T, build_vignettes=F)


library(tictoc)
library(FBMS)
#library(devtools)
#devtools::install_github("jonlachmann/irls.sgd", force=T, build_vignettes=F)
library(irls.sgd)
#Kaggle API 
library(RKaggle)

# Download latest version
df <- RKaggle::get_dataset("alexteboul/heart-disease-health-indicators-dataset")

summary(df)
dim(df)



#number of observations and covariates in the data

n <- dim(df)[1] 
p <- dim(df)[2] - 1   

params <- gen.params.gmjmcmc(p)
transforms <- c("sigmoid","pm1","p0","p05","p2","p3")

r = 0.01   # Parameter for the model prior

logistic.posterior.bic.irlssgd <- function (y, x, model, complex, mlpost_params) 
{
  if (!is.null(mlpost_params$crit)) {
    mod <- glm.sgd(x[,model], y, binomial(), 
                   sgd.ctrl = list(start=mlpost_params$coefs, subs=mlpost_params$subs, 
                                   maxit=10, alpha=0.00008, decay=0.99, histfreq=10))
    mod$deviance <- get_deviance(mod$coefficients, x[,model], y, binomial())
    mod$rank <- length(mod$coefficients)
  } else {
    mod <- irls.sgd(as.matrix(x[,model]), y, binomial(),
                    irls.control=list(subs=mlpost_params$subs, maxit=20, tol=1e-7, 
                                      cooling = c(1,0.9,0.75), expl = c(3,1.5,1)),
                    sgd.control=list(subs=mlpost_params$subs, maxit=250, alpha=0.001, 
                                     decay=0.99, histfreq=10))
  }
  
  # logarithm of marginal likelihood
  mloglik <- -mod$deviance / 2 - 0.5 * log(length(y)) * (mod$rank - 1)
  
  # logarithm of model prior
  if (length(mlpost_params$r) == 0)  mlpost_params$r <- 1/dim(x)[1]  # default value or parameter r
  lp <- log_prior(mlpost_params, complex)
  crit <- mloglik + lp
  
  if (!is.null(mlpost_params$crit) && mlpost_params$crit > crit) {
    return(list(crit = mlpost_params$crit, coefs = mlpost_params$coefs))
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
result1 <- fbms(formula = HeartDiseaseorAttack ~ 1 + ., data = df, P = 2,
                transforms = transforms, params = params, method = "gmjmcmc",
                family = "custom", loglik.pi = logistic.posterior.bic.irlssgd,  
                model_prior = list(r = r, subs = 0.01),  sub  = T)            
time1 <- toc() 

set.seed(100002)
# regular analysis
tic()
result2 <- fbms(formula = HeartDiseaseorAttack ~ 1 + ., data = df, P = 2,
                transforms = transforms, params = params, method = "gmjmcmc", 
                family = "binomial", beta_prior = list(type = "Jeffreys-BIC"),
                model_prior = list(r = r))
time2 <- toc() 

c(time1, time2)

############################
#
# More serious analysis
#
############################


# with subsampling

set.seed(100003)

tic()
result_parallel_1 <- fbms(formula = HeartDiseaseorAttack ~ 1 + ., data = df, P = 3,
                          transforms = transforms, params = params, 
                          method = "gmjmcmc.parallel", runs = 10, cores = 10,
                          family = "custom", loglik.pi = logistic.posterior.bic.irlssgd, 
                          model_prior = list(r = r, subs = 0.01), sub  = T)
time3 <- toc() 

summary(result_parallel_1)

# without subsampling


set.seed(100004)

tic()
result_parallel_2 <- fbms(formula = HeartDiseaseorAttack ~ 1 + ., data = df, P = 3, 
                          transforms = transforms, params = params, family = "binomial", 
                          method = "gmjmcmc.parallel", runs = 10, cores = 10,
                          model_prior = list(r = r),
                          beta_prior = list(type = "Jeffreys-BIC"))
time4 <- toc() 

summary(result_parallel_2)

filename = paste0("Ex11_Results_",r,"_4.RData")
save.image(filename)

############################
#
# Final analysis
#
############################

# with subsampling

set.seed(100005)
tic()
result_parallel_long_1 <- fbms(formula = HeartDiseaseorAttack ~ 1 + ., data = df, P = 10,
                          transforms = transforms, params = params, N = 500,
                          method = "gmjmcmc.parallel", runs = 40, cores = 40,
                          family = "custom", loglik.pi = logistic.posterior.bic.irlssgd, 
                          model_prior = list(r = r, subs = 0.01), sub  = T)

time5 <- toc() 
summary(result_parallel_long_1)

filename = paste0("Ex11_Results_",r,"_5.RData")
save.image(filename)

# regular analysis


set.seed(100006)

tic()
result_parallel_long_2 <- fbms(formula = HeartDiseaseorAttack ~ 1 + ., data = df, P = 10, 
                          transforms = transforms, params = params, family = "binomial", 
                          method = "gmjmcmc.parallel", runs = 40, cores = 40, N = 500,
                          model_prior = list(r = r),
                          beta_prior = list(type = "Jeffreys-BIC"))
time6 <- toc() 


summary(result_parallel_long_2)


############################################################################

C = cor(df, use = "everything",
        method = "spearman")

corrplot::corrplot(C)

apply((abs(C - diag(diag(C)))), 2, max)

filename = paste0("Ex11_Results_",r,".RData")
save.image(filename)