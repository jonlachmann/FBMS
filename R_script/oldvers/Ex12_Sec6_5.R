#######################################################
#
# Example 13 (Section 6.5):
#
# Cox Regression (using only fbms)
#
# This is the valid version for the JSS Paper
#
#######################################################

#install.packages("FBMS")
library(FBMS)
library(pec) #for the computation of cindex

#install.packages("survival")
library(survival)


# Download data
download.file('https://www.uniklinik-freiburg.de/fileadmin/mediapool/08_institute/biometrie-statistik/Dateien/Studium_und_Lehre/Lehrbuecher/Multivariable_Model-building/gbsg_br_ca.zip',
              'gbsg_br_ca.zip')
df1 <- read.csv(unz('gbsg_br_ca.zip',
                    'gbsg_br_ca/gbsg_br_ca.csv'),
                header = TRUE)
file.remove("gbsg_br_ca.zip")

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
probs$gen <- c(1,1,0,0)


# specify the custom function to estimate log posterior for cox 
surv.pseudo.loglik = function(y, x, model, complex, mlpost_params){
  
  data <- data.frame(time = mlpost_params$time, cens = y, as.matrix(x[,model]))[,-3]  # Removing intercept
  if(dim(data)[2]==2)
  {  
    return(list(crit=-.Machine$double.xmax, coefs=rep(0,1)))
  }  else {
    formula1 <- as.formula(paste0("Surv(time,cens)","~ 1 + ."))
    
    out <- coxph(formula1, data = data)
    
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
}



#######################################################
#
#   Analysis with 4 different modeling strategies
#
#


# 1) Single chain analysis (just to illustrate how it works)
set.seed(121)
result1 <- fbms(formula = cens ~ 1 + .,data = df.train[,-1], params = params, P = 5,
                transforms = transforms, method = "gmjmcmc",
                family = "custom", loglik.pi = surv.pseudo.loglik, 
                model_prior = list(r = 0.5),
                extra_params = list(time = time))


summary(result1,labels = names(df.train)[-(1:2)])

# 2) Parallel version only linear terms


set.seed(122)
result2 <- fbms(formula = cens ~ 1 + .,data = df.train[,-1], params = params,
                method = "mjmcmc.parallel", runs = 40, cores = 40,
                family = "custom", loglik.pi = surv.pseudo.loglik, 
                model_prior = list(r = 0.5), extra_params = list(time = time))

summary(result2,tol = 0.01,labels = names(df.train)[-(1:2)],effects = c(0.025,0.5,0.975))



# 3) Parallel version only fractional polynomials

set.seed(123)
probs$gen <- c(0,1,0,1)
params$feat$D <- 1

result3 <- fbms(formula = cens ~ 1 + .,data = df.train[,-1], params = params, probs = probs, P = 10,
                transforms = transforms, method = "gmjmcmc.parallel", runs = 40, cores = 40,
                family = "custom", loglik.pi = surv.pseudo.loglik, 
                model_prior = list(r = 0.5), extra_params = list(time = time))


summary(result3,tol = 0.01, effects = c(0.025,0.5,0.975))



# 4) Parallel version using all types of non-linear features
set.seed(124)
probs$gen <- c(1,1,1,1)
params$feat$D <- 5
result4 <- fbms(formula = cens ~ 1 + .,data = df.train[,-1], params = params, probs = probs,P = 20, 
                transforms = transforms, method = "gmjmcmc.parallel", runs = 40, cores = 40,
                family = "custom", loglik.pi = surv.pseudo.loglik, 
                model_prior = list(r = 0.5), extra_params = list(time = time))


summary(result4,tol = 0.01)






################################################
#
#  Prediction and C index using model averaging
#
################################################


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

# Clean the train and test data for the next type of predictions

df.train <- df[train,]
df.test <- df[-train,]

##############################################
#
# Prediction and C index using best model
#
##############################################


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

all.cindices <- rbind(all.cindices, round(unlist(c(cindex1, cindex2, cindex3, cindex4, cindex5, cindex6)),3))

# Clean the train and test data for the next type of predictions

df.train <- df[train,]
df.test <- df[-train,]

##############################################
#
# Prediction and C index using mpm model
#
##############################################


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


all.cindices <- rbind(all.cindices, round(unlist(c(cindex1, cindex2, cindex3, cindex4, cindex5, cindex6)),3))
rownames(all.cindices) = c("Model Averaging", "Best Model", "MPM")

print(all.cindices)
