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

library(FBMS)

n = 2000
p = 50

set.seed(1)
X2 <- as.data.frame(array(data = rbinom(n = n*p,size = 1,prob = runif(n = n*p,0,1)),dim = c(n,p)))
y2.Mean = 1+7*(X2$V4*X2$V17*X2$V30*X2$V10) + 9*(X2$V7*X2$V20*X2$V12)+ 3.5*(X2$V9*X2$V2)+1.5*(X2$V37)
Y2 <- rnorm(n = n,mean = y2.Mean,sd = 1)
df <- data.frame(Y2,X2)
summary(df)

str(df)


# Split data into training and test dataset
df.training <- df[1:(n/2),]
df.test <- df[(n/2 + 1):n,]
df.test$Mean <- y2.Mean[(n/2 + 1):n]



#############################################################################
#
#   FBMS logic regression with a Jeffreys parameter prior
#
#############################################################################



# FBMS - unlike the EMJMCMC package - does not explicitly have GMJMCMC for logic regression, 
# but we can easily run it without an "or" operator as "and" and "not" allow 
# to compute "or" through de Morgan law.

transforms <- c("not")
probs <- gen.probs.gmjmcmc(transforms)
probs$gen <- c(1,1,0,1) #No projections allowed

params <- gen.params.gmjmcmc(p)
params$feat$pop.max <- 50
params$feat$L <- 15



estimate.logic.lm = function(y, x, model, complex, mlpost_params)
{
  # Computation of marginal log-likelihood using Jeffreys prior
  suppressWarnings({
    mod <- fastglm(as.matrix(x[, model]), y, family = gaussian())
  })
  mloglik <- -(mod$aic + (log(length(y))-2) * (mod$rank))/2 
  
  # Computation of log of model prior
  wj <- complex$width
  lp <- sum(log(factorial(wj))) - sum(wj*log(4*mlpost_params$p) - log(4))
  
  # log posterior up to a constant
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

result <- fbms(formula = Y2~1+., data = df.training, probs = probs, params = params,  
               method = "gmjmcmc", transforms = transforms, N = 500, P = 25,
               family = "custom", loglik.pi = estimate.logic.lm,
               model_prior = list(p = p))
summary(result)
mpm <- get.mpm.model(result, y = df.training$Y2, x = df.training[,-1], family = "custom", loglik.pi = estimate.logic.lm,params = list(p = 50))
mpm$coefs
mpm <- get.mpm.model(result, y = df.training$Y2, x = df.training[,-1])
mpm$coefs
mbest <- get.best.model(result)
mbest$coefs


pred <- predict(result, x =  df.test[,-1], link = function(x)(x))  
pred_mpm <- predict(mpm, x =  df.test[,-1], link = function(x)(x))
pred_best <- predict(mbest, x =  df.test[,-1], link = function(x)(x))


#prediction errors
sqrt(mean((pred$aggr$mean - df.test$Y2)^2))
sqrt(mean((pred_mpm - df.test$Y2)^2))
sqrt(mean((pred_best - df.test$Y2)^2))
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

result_parallel <- fbms(formula = Y2~1+.,data = df.training, probs = probs, params = params, 
                   method = "gmjmcmc.parallel", transforms = transforms, N = 500, P=25,
                   family = "custom", loglik.pi = estimate.logic.lm, 
                   model_prior = list(p = p), runs = 16, cores = 8)
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


library(BAS) #needed for hypergeometric functions
estimate.logic.tcch = function(y, x, model, complex, mlpost_params)
{
  # Computation of marginal log likelihood
  
  suppressWarnings({
    mod <- fastglm(as.matrix(x[, model]), y, family = gaussian())
  })
  
  p.v <- (mlpost_params$n+1)/(mod$rank+1)
  
  y_mean <- mean(y)
  TSS <- sum((y - y_mean)^2)
  RSS <- sum(mod$residuals^2)
  R.2 <- 1 - (RSS / TSS)
  p <- mod$rank
  
  mloglik = (-0.5*p*log(p.v) -0.5*(mlpost_params$n-1)*log(1-(1-1/p.v)*R.2) + log(beta((mlpost_params$p.a+p)/2,mlpost_params$p.b/2)) - log(beta(mlpost_params$p.a/2,mlpost_params$p.b/2)) + log(phi1(mlpost_params$p.b/2,(mlpost_params$n-1)/2,(mlpost_params$p.a+mlpost_params$p.b+p)/2,mlpost_params$p.s/2/p.v,R.2/(p.v-(p.v-1)*R.2))) - hypergeometric1F1(mlpost_params$p.b/2,(mlpost_params$p.a+mlpost_params$p.b)/2,mlpost_params$p.s/2/p.v,log = T)) 
  if(mloglik ==-Inf||is.na(mloglik )||is.nan(mloglik ))
    mloglik  = -10000
 
  # Computation of log of model prior
  
  wj <- complex$width
  lp <- sum(log(factorial(wj))) - sum(wj*log(mlpost_params$p) + (2*wj-2)*log(2))
  
   
  logpost <- mloglik + lp + mlpost_params$n
  
  if(logpost==-Inf)
    logpost = -10000
  
  return(list(crit = logpost, coefs = mod$coefficients))
}


set.seed(5001)



result.tcch <- fbms(formula = Y2~1+.,data = df.training, probs = probs, params = params,
               method = "gmjmcmc", transforms = transforms, N = 500, P = 25,
               family = "custom", loglik.pi = estimate.logic.tcch,
               model_prior = list(p = p, n = n),
               beta_prior =  list(p.a = 1, p.b = 1, p.r = 1.5, p.s = 0, p.k = 1))
summary(result.tcch)
mpm.tcch <- get.mpm.model(result.tcch, y = df.training$Y2, x = df.training[,-1], family = "custom", loglik.pi = estimate.logic.lm,params = list(p = 50, n = n, p.a = 1, p.b = 1, p.r = 1.5, p.s = 0, p.k = 1))
mbest.tcch <- get.best.model(result.tcch)


pred.tcch <- predict(result.tcch, x =  df.test[,-1], link = function(x)(x))  
pred_mpm.tcch <- predict(mpm.tcch, x =  df.test[,-1], link = function(x)(x))
pred_best.tcch <- predict(mbest.tcch, x =  df.test[,-1], link = function(x)(x))


#prediction errors
sqrt(mean((pred.tcch$aggr$mean - df.test$Y2)^2))
sqrt(mean((pred_best.tcch - df.test$Y2)^2))
sqrt(mean((pred_mpm.tcch - df.test$Y2)^2))
sqrt(mean((df.test$Mean - df.test$Y2)^2))

#prediction errors to the true means
sqrt(mean((pred.tcch$aggr$mean - df.test$Mean)^2))
sqrt(mean((pred_best.tcch - df.test$Mean)^2))
sqrt(mean((pred_mpm.tcch - df.test$Mean)^2))



plot(pred.tcch$aggr$mean, df.test$Y2)
points(pred.tcch$aggr$mean,df.test$Mean,col = 2)
points(pred_best.tcch,df.test$Mean,col = 3)
points(pred_mpm.tcch,df.test$Mean,col = 4)


# Now parallel inference

set.seed(5002)

result_parallel.tcch <- fbms(formula = Y2~1+.,data = df.training, probs = probs, params = params,
                    method = "gmjmcmc.parallel", transforms = transforms, N = 500, P = 25,
                    family = "custom", loglik.pi = estimate.logic.tcch,
                    runs = 16, cores = 8, model_prior = list(p = p, n = n),
                    beta_prior =  list(p.a = 1, p.b = 1, p.r = 1.5, p.s = 0, p.k = 1))
summary(result_parallel.tcch)
mpm <- get.mpm.model(result_parallel.tcch,y = df.training$Y2,x = df.training[,-1],family = "custom", loglik.pi = estimate.logic.lm,params = list(p = 50, n = n, p.a = 1, p.b = 1, p.r = 1.5, p.s = 0, p.k = 1))
mbest <- get.best.model(result_parallel.tcch)


pred_parallel.tcch <- predict(result_parallel.tcch, x =  df.test[,-1], link = function(x)(x))  
pred_par_mpm.tcch <- predict(mpm, x =  df.test[,-1], link = function(x)(x))
pred_par_best.tcch <- predict(mbest, x =  df.test[,-1], link = function(x)(x))


#prediction errors
sqrt(mean((pred_parallel.tcch$aggr$mean - df.test$Y2)^2))
sqrt(mean((pred_par_best.tcch - df.test$Y2)^2))
sqrt(mean((pred_par_mpm.tcch - df.test$Y2)^2))
sqrt(mean((df.test$Mean - df.test$Y2)^2))

#prediction errors to the true means
sqrt(mean((pred_parallel.tcch$aggr$mean - df.test$Mean)^2))
sqrt(mean((pred_par_best.tcch - df.test$Mean)^2))
sqrt(mean((pred_par_mpm.tcch - df.test$Mean)^2))



plot(pred_parallel.tcch$aggr$mean, df.test$Y2)
points(pred_parallel.tcch$aggr$mean,df.test$Mean,col = 2)
points(pred_par_best.tcch,df.test$Mean,col = 3)
points(pred_par_mpm.tcch,df.test$Mean,col = 4)