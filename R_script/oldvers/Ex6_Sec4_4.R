#######################################################
#
# Example 6 (Section 4.4):
#
# Prediction using non-linear Projections, only using fbms function
#
#  DATA - abalone data set
#
#  Data set is available at https://www.kaggle.com/datasets/rodolfomendes/abalone-dataset
#
#  For convenience we provide the file abalone.csv which contains already the names of the variables
#
#
# This is the valid version for the JSS Paper
#
#######################################################


library(FBMS)

data("abalone")

df = abalone
df$Sex_F_vs_I = as.numeric(df$Sex == "F")
df$Sex_M_vs_I = as.numeric(df$Sex == "M")
df$Sex = as.factor(df$Sex)
df$Rings = as.numeric(df$Rings)

summary(df)


#number of observations in the data

n = dim(df)[1] 

# Create Training and Test dataset

# Sam Waugh (1995) "Extending and benchmarking Cascade-Correlation", PhD
# thesis, Computer Science Department, University of Tasmania.

#-- Test set performance (final 1044 examples, first 3133 used for training):
#  24.86% Cascade-Correlation (no hidden nodes)
#  26.25% Cascade-Correlation (5 hidden nodes)
#  21.5%  C4.5
#  0.0%  Linear Discriminate Analysis
#  3.57% k=5 Nearest Neighbour
#  (Problem encoded as a classification task)


# remove variable sex because gmjmcmc cannot handle factor variables
df.training = df[1:3133,-2]
df.test = df[3134:n,-2]


summary(df.training)


pred.RMSE = rep(0,5)   # To collect the results of prediction RMSE from the five different methods

pred.RMSE.mpm = rep(0,5) # Same for MPM
pred.RMSE.best = rep(0,5) # Same for posterior modes in the model space

transforms = c("sigmoid")
probs = gen.probs.gmjmcmc(transforms)
probs$gen = c(0,0,1,1) #Only projections!


#############################################################################
#
#   Using method 0 for alpha (simply set to 1, default)
#
#############################################################################

set.seed(5001)
result = fbms(data = df.training, method = "gmjmcmc", transforms = transforms, 
                probs = probs)
summary(result)

pred = predict(result, x =  df.test[,-1])  
pred.RMSE[1] = sqrt(mean((pred$aggr$mean - df.test$Rings)^2))

preds = predict(get.best.model(result), df.test[, -1])
pred.RMSE.best[1] = sqrt(mean((preds - df.test$Rings)^2))

preds = predict(get.mpm.model(result, y = df.training$Rings, x = df.training[, -1]), df.test[, -1])
pred.RMSE.mpm[1] = sqrt(mean((preds - df.test$Rings)^2))


plot(pred$aggr$mean, df.test$Rings)





#############################################################################
#
#   Parallel version
#
#############################################################################

set.seed(5003)
result_parallel = fbms(data = df.training, method = "gmjmcmc.parallel", runs = 40, cores = 40, 
                         transforms = transforms, probs = probs, P=25)
summary(result_parallel, tol = 0.05)



pred_parallel = predict(result_parallel, x =  df.test[,-1], link = function(x)(x))  
pred.RMSE[2] = sqrt(mean((pred_parallel$aggr$mean - df.test$Rings)^2))

preds = predict(get.best.model(result_parallel), df.test[, -1])
pred.RMSE.best[2] = sqrt(mean((preds - df.test$Rings)^2))

preds = predict(get.mpm.model(result_parallel, y = df.training$Rings, x = df.training[, -1]), df.test[, -1])
pred.RMSE.mpm[2] = sqrt(mean((preds - df.test$Rings)^2))



plot(pred_parallel$aggr$mean, df.test$Rings)
abline(0,1)



#############################################################################
#
#   Using method 3 to estimate alpha
#
#############################################################################


params = gen.params.gmjmcmc(ncol(df.training) - 1)
params$feat$alpha = "deep"


set.seed(5003)
result.a3 = fbms(data = df.training, method = "gmjmcmc", transforms = transforms, 
                   probs = probs, params = params)
summary(result.a3)



pred.a3 = predict(result.a3, x =  df.test[,-1], link = function(x)(x))  
pred.RMSE[3] = sqrt(mean((pred.a3$aggr$mean - df.test$Rings)^2))

plot(pred.a3$aggr$mean, df.test$Rings)


preds = predict(get.best.model(result.a3), df.test[, -1])
pred.RMSE.best[3] = sqrt(mean((preds - df.test$Rings)^2))

#not yet applicable to the deep method
#preds = predict(get.mpm.model(result.a3, y = df.training$Rings, x = df.training[, -1]), df.test[, -1])
pred.RMSE.mpm[3] =  NA#sqrt(mean((preds - df.test$Rings)^2))




#############################################################################
#
#   Parallel version  params$feat$alpha = "random"
#
#############################################################################

params$feat$alpha = "random"

set.seed(5004)
result_parallel.a3 = fbms(data = df.training, method = "gmjmcmc.parallel",  runs = 40, cores = 40,
                            transforms = transforms, probs = probs, params = params, P=25)
summary(result_parallel.a3, tol = 0.05)



pred_parallel.a3 = predict(result_parallel.a3, x =  df.test[,-1], link = function(x)(x))  
pred.RMSE[4] = sqrt(mean((pred_parallel.a3$aggr$mean - df.test$Rings)^2))

preds = predict(get.best.model(result_parallel.a3), df.test[, -1])
pred.RMSE.best[4] = sqrt(mean((preds - df.test$Rings)^2))

#preds = predict(get.mpm.model(result_parallel, y = df.training$Rings, x = df.training[, -1]), df.test[, -1])
pred.RMSE.mpm[4] = NA #sqrt(mean((preds - df.test$Rings)^2))


plot(pred_parallel.a3$aggr$mean, df.test$Rings)
abline(0,1)



#############################################################################
#
#   Parallel version with fractional polynomials
#
#############################################################################

transforms = c("p0","p2","p3","p05","pm05","pm1","pm2","p0p0","p0p05","p0p1","p0p2","p0p3","p0p05","p0pm05","p0pm1","p0pm2")
probs = gen.probs.gmjmcmc(transforms)
probs$gen = c(0,1,0,1) #Only modifications!

set.seed(50005)
result.fp = fbms(data = df.training, method = "gmjmcmc.parallel", runs = 40, cores = 40, 
                   transforms = transforms, probs = probs, P=25)
summary(result.fp)


pred_fp = predict(result.fp, x =  df.test[,-1], link = function(x)(x))  
pred.RMSE[5] = sqrt(mean((pred_fp$aggr$mean - df.test$Rings)^2))

plot(pred_fp$aggr$mean, df.test$Rings)

preds = predict(get.best.model(result.fp), df.test[, -1])
pred.RMSE.best[5] = sqrt(mean((preds - df.test$Rings)^2))

preds = predict(get.mpm.model(result.fp, y = df.training$Rings, x = df.training[, -1]), df.test[, -1])
pred.RMSE.mpm[5] = sqrt(mean((preds - df.test$Rings)^2))



#############################################################################
#
#   Summary of predictions
#
#############################################################################

round(pred.RMSE,3)
round(pred.RMSE.best,3)
round(pred.RMSE.mpm,3)

