#######################################################
#
# Example 7 (Section 5.2):
#
# Logic regression with a different model prior
#
# DATA - simulated 
#
#
#
# This is the valid version for the JSS Paper
#
#######################################################

#library(devtools)
#devtools::install_github("jonlachmann/GMJMCMC@FBMS", force=T, build_vignettes=F)

library(FBMS)
use.fbms <- TRUE  

set.seed(1)
X2 <- as.data.frame(array(data = rbinom(n = 50*2000,size = 1,prob = runif(n = 50*2000,0,1)),dim = c(2000,50)))
Y2 <- rnorm(n = 2000,mean = 1+7*(X2$V4*X2$V17*X2$V30*X2$V10) + 9*(X2$V7*X2$V20*X2$V12)+ 3.5*(X2$V9*X2$V2)+1.5*(X2$V37),sd = 1)
df <- data.frame(Y2,X2)
summary(df)

str(df)

#number of observations in the data

n = dim(df)[1] 



# remove variable sex because gmjmcmc cannot handle factor variables
df.training <- df[1:1000,]
df.test <- df[1001:n,]
df.test$Mean <- (1+7*(X2$V4*X2$V17*X2$V30*X2$V10) + 9*(X2$V7*X2$V20*X2$V12)+ 3.5*(X2$V9*X2$V2)+1.5*(X2$V37))[1001:n]

#FBMS unlike EMJMCMC package does not explicitly have GMJMCMC for logic regression, but we can easily run it without an 
# "or" operator as "and" and "not" allow to automatically handle "or" through de Morgan law.
transforms <- c("not")
probs <- gen.probs.gmjmcmc(transforms)
probs$gen <- c(1,1,0,1) #No projections allowed


params <- gen.params.gmjmcmc(ncol(df.training) - 1)
params$feat$pop.max <- 51
params$feat$L <- 15
##############################################
#############################################################################
#
#   FBMS logic regression with a Jeffreys parameter prior
#
#############################################################################


estimate.logic.lm = function(y, x, model, complex, mlpost_params)
{
  
  if (length(mlpost_params) == 0) 
    mlpost_params <- list(r = 1/dim(x)[1]) 
  suppressWarnings({
    mod <- fastglm(as.matrix(x[, model]), y, family = gaussian())
  })
  
  wj <- complex$width
 
  lp <- sum(log(factorial(wj))) - sum(wj*log(mlpost_params$p) + (2*wj-2)*log(2))
  
  mloglik <- -(mod$aic + (log(length(y))-2) * (mod$rank))/2 
  
  logpost <- mloglik + lp 
  
  if(logpost==-Inf)
    logpost = -10000
  
  return(list(crit = logpost, coefs = mod$coefficients))
}



#############################################################################
#
#   Logic regression training
#
#############################################################################

set.seed(5001)

if (use.fbms) {
  result <- fbms(formula = Y2~1+., data = df.training, family = "custom", loglik.pi = estimate.logic.lm,N = 500,N.final = 500, P = 25,
                 method = "gmjmcmc", model_prior = list(p = 50), beta_prior = NULL, transforms = transforms, 
                 probs = probs, params = params)
} else {
  #  result <- gmjmcmc(df.training, transforms = transforms, probs = probs)
  
  result <- gmjmcmc(x = df.training[, -1], y = df.training[, 1], loglik.pi = estimate.logic.lm,N = 500,N.final = 500, , P = 25,
                    transforms = transforms, mlpost_params = list(p = 50), params = params, probs = probs)
  
}
summary(result)
mpm <- get.mpm.model(result, y = df.training$Y2, x = df.training[,-1], family = "custom", loglik.pi = estimate.logic.lm,params = list(p = 50))
mbest <- get.best.model(result)


pred <- predict(result, x =  df.test[,-1], link = function(x)(x))  
pred_mpm <- predict(mpm, x =  df.test[,-1], link = function(x)(x))
pred_best <- predict(mbest, x =  df.test[,-1], link = function(x)(x))


#prediction errors
sqrt(mean((pred$aggr$mean - df.test$Y2)^2))
sqrt(mean((pred_best - df.test$Y2)^2))
sqrt(mean((pred_mpm - df.test$Y2)^2))
sqrt(mean((df.test$Mean - df.test$Y2)^2))

#prediction errors to the true means
sqrt(mean((pred$aggr$mean - df.test$Mean)^2))
sqrt(mean((pred_best - df.test$Mean)^2))
sqrt(mean((pred_mpm - df.test$Mean)^2))



plot(pred$aggr$mean, df.test$Y2)
points(pred$aggr$mean,df.test$Mean,col = 2)
points(pred_best,df.test$Mean,col = 3)
points(pred_mpm,df.test$Mean,col = 4)




#############################################################################
#
#   Parallel version just 16 chains on 8 cores
#
#############################################################################


set.seed(5002)

if (use.fbms) {
  result_parallel <- fbms(formula = Y2~1+.,data = df.training, family = "custom", loglik.pi = estimate.logic.lm, N = 500, N.final = 500,
                          method = "gmjmcmc.parallel",model_prior = list(p = 50), beta_prior = NULL, runs = 16, cores = 8,
                          transforms = transforms, probs = probs, params = params, P=25)
} else {
  result_parallel =  gmjmcmc.parallel(runs = 16, cores = 8, x = df.training[, -1], y = df.training[, 1],
                                      loglik.pi = estimate.logic.lm, mlpost_params = list(p = 50), N = 500,N.final = 500,
                                      transforms = transforms, probs = probs, params = params, P=25)
}
summary(result_parallel)
mpm <- get.mpm.model(result_parallel, y = df.training$Y2, x = df.training[,-1], family = "custom", loglik.pi = estimate.logic.lm,params = list(p = 50))
mbest <- get.best.model(result_parallel)


pred_parallel <- predict(result_parallel, x =  df.test[,-1], link = function(x)(x))  
pred_par_mpm <- predict(mpm, x =  df.test[,-1], link = function(x)(x))
pred_par_best <- predict(mbest, x =  df.test[,-1], link = function(x)(x))


#prediction errors
sqrt(mean((pred_parallel$aggr$mean - df.test$Y2)^2))
sqrt(mean((pred_par_best - df.test$Y2)^2))
sqrt(mean((pred_par_mpm - df.test$Y2)^2))
sqrt(mean((df.test$Mean - df.test$Y2)^2))

#prediction errors to the true means
sqrt(mean((pred_parallel$aggr$mean - df.test$Mean)^2))
sqrt(mean((pred_par_best - df.test$Mean)^2))
sqrt(mean((pred_par_mpm - df.test$Mean)^2))



plot(pred_parallel$aggr$mean, df.test$Y2)
points(pred_parallel$aggr$mean,df.test$Mean,col = 2)
points(pred_par_best,df.test$Mean,col = 3)
points(pred_par_mpm,df.test$Mean,col = 4)



#############################################################################
#
#   FBMS logic regression with a tCCH parameter prior
#
#############################################################################

transforms <- c("not")
probs <- gen.probs.gmjmcmc(transforms)
probs$gen <- c(1,1,0,1) #No projections allowed
probs$filter <- 0.6
params <- gen.params.gmjmcmc(ncol(df.training) - 1)
params$feat$pop.max <- 51

library(BAS) #needed for hypergeometric functions
estimate.logic.tcch = function(y, x, model, complex, mlpost_params)
{
  
  if (length(mlpost_params) == 0) 
    mlpost_params <- list(r = 1 / dim(x)[1])
  suppressWarnings({
    mod <- fastglm(as.matrix(x[, model]), y, family = gaussian())
  })
  
  wj <- complex$width
  
  lp <- sum(log(factorial(wj))) - sum(wj*log(mlpost_params$p) + (2*wj-2)*log(2))
  
  p.v <- (mlpost_params$n+1)/(mod$rank+1)
  
  y_mean <- mean(y)
  TSS <- sum((y - y_mean)^2)
  RSS <- sum(mod$residuals^2)
  R.2 <- 1 - (RSS / TSS)
  p <- mod$rank
  
  mloglik = (-0.5*p*log(p.v) -0.5*(mlpost_params$n-1)*log(1-(1-1/p.v)*R.2) + log(beta((mlpost_params$p.a+p)/2,mlpost_params$p.b/2)) - log(beta(mlpost_params$p.a/2,mlpost_params$p.b/2)) + log(phi1(mlpost_params$p.b/2,(mlpost_params$n-1)/2,(mlpost_params$p.a+mlpost_params$p.b+p)/2,mlpost_params$p.s/2/p.v,R.2/(p.v-(p.v-1)*R.2))) - hypergeometric1F1(mlpost_params$p.b/2,(mlpost_params$p.a+mlpost_params$p.b)/2,mlpost_params$p.s/2/p.v,log = T)) 
  if(mloglik ==-Inf||is.na(mloglik )||is.nan(mloglik ))
    mloglik  = -10000
  
  logpost <- mloglik + lp + mlpost_params$n
  
  if(logpost==-Inf)
    logpost = -10000
  
  return(list(crit = logpost + lp, coefs = mod$coefficients))
}


set.seed(5001)

if (use.fbms) {
  result <- fbms(formula = Y2~1+.,data = df.training, family = "custom", loglik.pi = estimate.logic.tcch,N = 500,N.final = 500, P = 25,
                 method = "gmjmcmc", transforms = transforms, 
                 probs = probs,model_prior = list(p = 50,n = n),beta_prior =  list(p.a = 1, p.b = 1, p.r = 1.5, p.s = 0, p.k = 1), params = params)
} else {
  #  result <- gmjmcmc(df.training, transforms = transforms, probs = probs)
  
  result <- gmjmcmc(x = df.training[, -1], y = df.training[, 1], loglik.pi = estimate.logic.tcch,N = 500,N.final = 500, P = 25,
                    transforms = transforms,mlpost_params = list(p = 50, n = n, p.a = 1, p.b = 1, p.r = 1.5, p.s = 0, p.k = 1), params = params, probs = probs)
  
}
summary(result)
mpm <- get.mpm.model(result, y = df.training$Y2, x = df.training[,-1], family = "custom", loglik.pi = estimate.logic.lm,params = list(p = 50, n = n, p.a = 1, p.b = 1, p.r = 1.5, p.s = 0, p.k = 1))
mbest <- get.best.model(result)


pred <- predict(result, x =  df.test[,-1], link = function(x)(x))  
pred_mpm <- predict(mpm, x =  df.test[,-1], link = function(x)(x))
pred_best <- predict(mbest, x =  df.test[,-1], link = function(x)(x))


#prediction errors
sqrt(mean((pred$aggr$mean - df.test$Y2)^2))
sqrt(mean((pred_best - df.test$Y2)^2))
sqrt(mean((pred_mpm - df.test$Y2)^2))
sqrt(mean((df.test$Mean - df.test$Y2)^2))

#prediction errors to the true means
sqrt(mean((pred$aggr$mean - df.test$Mean)^2))
sqrt(mean((pred_best - df.test$Mean)^2))
sqrt(mean((pred_mpm - df.test$Mean)^2))



plot(pred$aggr$mean, df.test$Y2)
points(pred$aggr$mean,df.test$Mean,col = 2)
points(pred_best,df.test$Mean,col = 3)
points(pred_mpm,df.test$Mean,col = 4)


# Now parallel inference

set.seed(5002)

if (use.fbms) {
  result_parallel <- fbms(formula = Y2~1+.,data = df.training, family = "custom", loglik.pi = estimate.logic.tcch,N = 500,N.final = 500,
                          method = "gmjmcmc.parallel", runs = 16, cores = 8, model_prior = list(p = 50,n = n),beta_prior =  list(p.a = 1, p.b = 1, p.r = 1.5, p.s = 0, p.k = 1),
                          transforms = transforms, probs = probs, params = params, P=25)
} else {
  result_parallel =  gmjmcmc.parallel(runs = 16, cores = 8, x = df.training[, -1], y = df.training[, 1],
                                      loglik.pi = estimate.logic.tcch,N = 500,N.final = 500,
                                      mlpost_params = list(p = 50, n = n, p.a = 1, p.b = 1, p.r = 1.5, p.s = 0, p.k = 1),
                                      transforms = transforms, probs = probs, params = params, P=25)
}
summary(result_parallel)
mpm <- get.mpm.model(result_parallel,y = df.training$Y2,x = df.training[,-1],family = "custom", loglik.pi = estimate.logic.lm,params = list(p = 50, n = n, p.a = 1, p.b = 1, p.r = 1.5, p.s = 0, p.k = 1))
mbest <- get.best.model(result_parallel)


pred_parallel <- predict(result_parallel, x =  df.test[,-1], link = function(x)(x))  
pred_par_mpm <- predict(mpm, x =  df.test[,-1], link = function(x)(x))
pred_par_best <- predict(mbest, x =  df.test[,-1], link = function(x)(x))


#prediction errors
sqrt(mean((pred_parallel$aggr$mean - df.test$Y2)^2))
sqrt(mean((pred_par_best - df.test$Y2)^2))
sqrt(mean((pred_par_mpm - df.test$Y2)^2))
sqrt(mean((df.test$Mean - df.test$Y2)^2))

#prediction errors to the true means
sqrt(mean((pred_parallel$aggr$mean - df.test$Mean)^2))
sqrt(mean((pred_par_best - df.test$Mean)^2))
sqrt(mean((pred_par_mpm - df.test$Mean)^2))



plot(pred_parallel$aggr$mean, df.test$Y2)
points(pred_parallel$aggr$mean,df.test$Mean,col = 2)
points(pred_par_best,df.test$Mean,col = 3)
points(pred_par_mpm,df.test$Mean,col = 4)