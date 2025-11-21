#######################################################
#
# Example 8 (Section 6.1):
#
# Logistic Regression, using only fbms
#
# This is the valid version for the JSS Paper
#
#######################################################

library(FBMS)


library(kernlab)
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
#
# single thread analysis
#
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
#
# multiple thread analysis
#
####################################################

set.seed(6002)

result_parallel <- fbms(formula = y~1+.,data = df, method = "gmjmcmc.parallel", 
                        family = "binomial", beta_prior = list(type = "Jeffreys-BIC"),
                          runs = 16, cores = 16, transforms = transforms, 
                          probs = probs, params = params, P=25)

summary(result_parallel)

#######################
#
# Prediction accuracy
# IMPORTANT: specify correct link function for predict
#

# Model averaging
pred_parallel <- predict(result_parallel, x =  df[,-1], link = function(x)(1/(1+exp(-x))))  
mean(round(pred_parallel$aggr$mean)==df$y)

# Best Model
bm_parallel <- get.best.model(result_parallel)
pred_bm_parallel <-  predict(bm_parallel, df[,-1],link = function(x)(1/(1+exp(-x))))
mean(round(pred_bm_parallel)==df$y)

# Median Probability Model
mpm_parallel <-  predict(get.mpm.model(result = result_parallel,family = "binomial",y = df$y,x=df[,-1]), df[,-1],link = function(x)(1/(1+exp(-x))))
mean(round(mpm_parallel)==df$y)
