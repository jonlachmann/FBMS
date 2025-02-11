#######################################################
#
# Example 13 (Section 6.5):
#
# Cox Regression
#
# This is the valid version for the JSS Paper
#
#######################################################

library(devtools)
devtools::install_github("jonlachmann/GMJMCMC@FBMS", force=T, build_vignettes=F)

#install.packages("FBMS")
library(FBMS)
library(pec) #for the computation of cindex

#install.packages("survival")
library(survival)

setwd("/home/florian/FBMS/")

download.file('https://www.uniklinik-freiburg.de/fileadmin/mediapool/08_institute/biometrie-statistik/Dateien/Studium_und_Lehre/Lehrbuecher/Multivariable_Model-building/gbsg_br_ca.zip',
              'gbsg_br_ca.zip')
df1 <- read.csv(unz('gbsg_br_ca.zip',
                    'gbsg_br_ca/gbsg_br_ca.csv'),
                header = TRUE)
#system('rm whitehall1.zip')


df <- df1[, c(13, 14, 2:4, 6:8, 10:12)]
names(df) = c("time","cens",names(df)[3:ncol(df)])


set.seed(123)
train <- c(sample((1:nrow(df))[df$cens == 1], sum(df$cens)*2/3), # split separately events
           sample((1:nrow(df))[df$cens == 0], sum(!df$cens)*2/3)) # and censored observations


df.train <- df[train,]
df.test <- df[-train,]

time <- df.train$time


params <- gen.params.gmjmcmc(data = df.train[,-1])
params$loglik$r = 0.5 
params$loglik$time = time    #the time variable goes into the params structure

params$feat$keep.min = 0.2

transforms <- c("p0","p2","p3","p05","pm05","pm1","pm2")
probs <- gen.probs.gmjmcmc(transforms)
#probs$gen <- c(1,1,0,1)


#specify the estimator function for cox 
surv.pseudo.loglik = function(y, x, model, complex, params){
  
  if(length(params$r)==0)
    params$r = 0.5
  data <- data.frame(time = params$time, cens = y, as.matrix(x[,model]))[,-3]  # Removing intercept
  if(dim(data)[2]==2)
  {  
    return(list(crit=-10000, coefs=rep(0,1)))
  }  else {
     formula1 = as.formula(paste0("Surv(time,cens)","~ 1 + ."))

     out = coxph(formula1, data = data)

     # logarithm of marginal likelihood
     mloglik <- (out$loglik[2] - out$loglik[1]) -  log(length(y)) * (dim(data)[2] - 2)/2   
     
     # logarithm of model prior
     if (length(params$r) == 0)  params$r <- 1/dim(x)[1]  # default value or parameter r
     lp <- log.prior(params, complex)
     
     return(list(crit = mloglik + lp, coefs =  c(0,out$coefficients)))
     
  }
  
}

#Single chain analysis (just to illustrate that it works)
set.seed(5)
result <- gmjmcmc(data = df.train[,-1], loglik.pi = surv.pseudo.loglik, transforms = transforms, params = params,  P = 5)

summary(result,tol = 0.01,labels = names(df.train)[-c(1:2)],effects = c(0.025,0.5,0.975))
summary(result)
summary(result,tol = 0.01,labels = names(df.train)[-c(1:2)])

linpreds.train <- predict(result,df.train[,-(1:2)], link = function(x) x)
linpreds <- predict(result,df.test[,-(1:2)], link = function(x) x)

#plot(linpreds$aggr$mean)

#Parallel version
set.seed(15)
probs$gen <- c(1,1,1,1)
result2 <- gmjmcmc.parallel(runs = 80, cores = 40, data = df.train[,-1], 
                            loglik.pi = surv.pseudo.loglik, transforms = transforms, 
                            probs = probs, params = params,  P = 25)
summary(result2,tol = 0.01,labels = names(df.train)[-1],effects = c(0.025,0.5,0.975))

linpreds2.train <- predict(result2,df.train[,-(1:2)], link = function(x) x)
linpreds2 <- predict(result2,df.test[,-(1:2)], link = function(x) x)
#plot(linpreds2$aggr$mean)



#############################################
#Parallel version only linear terms
set.seed(25)
probs$gen <- c(0,0,0,1)
result3 <- gmjmcmc.parallel(runs = 80, cores = 40, data = df.train[,-1], 
                            loglik.pi = surv.pseudo.loglik, transforms = transforms, 
                            probs = probs, params = params,  P = 25)

summary(result3,tol = 0.01,labels = names(df.train)[-(1:2)],effects = c(0.025,0.5,0.975))

linpreds3.train <- predict(result3,df.train[,-(1:2)], link = function(x) x)
linpreds3 <- predict(result3,df.test[,-(1:2)], link = function(x) x)
#plot(linpreds2$aggr$mean)

#############################################
#Parallel version only fractional polynomials
set.seed(35)
probs$gen <- c(0,1,0,1)
result4 <- gmjmcmc.parallel(runs = 80, cores = 40, data = df.train[,-1], 
                            loglik.pi = surv.pseudo.loglik, transforms = transforms, 
                            probs = probs, params = params,  P = 25)

summary(result4,tol = 0.01,labels = names(df.train)[-(1:2)],effects = c(0.025,0.5,0.975))

linpreds4.train <- predict(result4,df.train[,-(1:2)], link = function(x) x)
linpreds4 <- predict(result4,df.test[,-(1:2)], link = function(x) x)
#plot(linpreds2$aggr$mean)

# Compute cindex using package pec

df.train$average.lin.pred1 <- linpreds.train$aggr$mean
df.train$average.lin.pred2 <- linpreds2.train$aggr$mean
df.train$average.lin.pred3 <- linpreds3.train$aggr$mean
df.train$average.lin.pred4 <- linpreds4.train$aggr$mean

df.test$average.lin.pred1 <- linpreds$aggr$mean
df.test$average.lin.pred2 <- linpreds2$aggr$mean
df.test$average.lin.pred3 <- linpreds3$aggr$mean
df.test$average.lin.pred4 <- linpreds4$aggr$mean

mod1 <- coxph(Surv(time, cens) ~ average.lin.pred1, data = as.data.frame(df.train), x = TRUE)
cindex1 <- cindex(mod1, mod1$formula, data = as.data.frame(df.test), cens.model = 'cox')$AppCindex

mod2 <- coxph(Surv(time, cens) ~ average.lin.pred2, data = as.data.frame(df.train), x = TRUE)
cindex2 <- cindex(mod2, mod2$formula, data = as.data.frame(df.test), cens.model = 'cox')$AppCindex

mod3 <- coxph(Surv(time, cens) ~ average.lin.pred3, data = as.data.frame(df.train), x = TRUE)
cindex3 <- cindex(mod3, mod3$formula, data = as.data.frame(df.test), cens.model = 'cox')$AppCindex

mod4 <- coxph(Surv(time, cens) ~ average.lin.pred4, data = as.data.frame(df.train),  x = TRUE)
cindex4 <- cindex(mod4, mod4$formula, data = as.data.frame(df.test), cens.model = 'cox')$AppCindex


#Full model
mod5 <- coxph(Surv(time, cens) ~ 1+., data = as.data.frame(df.train[,1:11]),x = T)
cindex5 <- cindex(mod5, mod5$formula, data = as.data.frame(df.test), cens.model = 'cox')$AppCindex

#Model without predictors
mod6 <- coxph(Surv(time, cens) ~ 1, data = as.data.frame(df.train[,1:11]),x = T)
cindex6 <- cindex(mod6, mod6$formula, data = as.data.frame(df.test), cens.model = 'cox')$AppCindex

round(unlist(c(cindex1, cindex2, cindex3, cindex4, cindex5, cindex6)),3)


