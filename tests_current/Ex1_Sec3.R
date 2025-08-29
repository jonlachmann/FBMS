#################################################
#
# Example 1:  
#
# Kepler Example with the most recent database update, only using fbms function
#
# This is the valid version for the JSS paper
#
##################################################


#install.packages("FBMS")
#install.packages("devtools")
#library(devtools)
#devtools::install_github("jonlachmann/GMJMCMC@data-inputs", force=T, build_vignettes=F)
library(FBMS)

data(exoplanet)

train.indx <- 1:500
df.train = exoplanet[train.indx, ]
df.test = exoplanet[-train.indx, ]


to3 <- function(x) x^3
transforms <- c("sigmoid","sin_deg","exp_dbl","p0","troot","to3")


####################################################
#
# single thread analysis (default values, Section 3.1)
#
####################################################


set.seed(123)

result.default <- fbms(formula = semimajoraxis ~ 1 + . , data = df.train, method = "gmjmcmc", transforms = transforms)



####################################################
#
# single thread analysis (more iterations, Section 3.2)
#
####################################################


set.seed(123)

result.P50 <- fbms(data = df.train, method = "gmjmcmc", transforms = transforms,
                     P = 50, N = 1000, N.final = 5000)

 
####################################################
#
# multiple thread analysis (Section 3.3)
#
####################################################

set.seed(123)

result_parallel <- fbms(data = df.train, method = "gmjmcmc.parallel", transforms = transforms,
                          runs = 40, cores = 40, P = 25)


####################################################
#
# Inspection of Results (Section 3.4)
#
####################################################

######################
# summary

summary(result.default)
summary(result.default, pop = "all", labels = paste0("x",1:length(df.train[,-1])))


summary(result.P50)
summary(result.P50, pop = "best", labels = paste0("x",1:length(df.train[,-1])))
summary(result.P50, pop = "last", labels = paste0("x",1:length(df.train[,-1])))
summary(result.P50, pop = "last", tol = 0.01, labels = paste0("x",1:length(df.train[,-1])))
summary(result.P50, pop = "all")

summary(result_parallel)
library(tictoc)
tic()
summary(result_parallel, tol = 0.01, pop = "all")
toc()




######################
# plot

pdf("result.pdf") 
plot(result.default)
dev.off()

plot(result.default)



pdf("result.P50.pdf") 
plot(result.P50)
dev.off()

plot(result.P50)



pdf("result_parallel.pdf") 
plot(result_parallel)
dev.off()

plot(result_parallel)


######################
# Prediction


#preds <- predict(result.default, df.test[,-1], link = function(x) x)
preds <-  predict(result.default, df.test[,-1])
rmse.default <- sqrt(mean((preds$aggr$mean - df.test$semimajoraxis)^2))

pdf("prediction.pdf") 
plot(preds$aggr$mean, df.test$semimajoraxis)
dev.off()

plot(preds$aggr$mean, df.test$semimajoraxis)






###############################


#preds.P50 = predict(result.P50, df.test[,-1], link = function(x) x)  
preds.P50 = predict(result.P50, df.test[,-1])  
rmse.P50 <-  sqrt(mean((preds.P50$aggr$mean - df.test$semimajoraxis)^2))

pdf("prediction.P50.pdf") 
plot(preds.P50$aggr$mean, df.test$semimajoraxis)
dev.off()

plot(preds.P50$aggr$mean, df.test$semimajoraxis)



###############################


preds.multi <- predict(result_parallel , df.test[,-1], link = function(x) x)
rmse.parallel <- sqrt(mean((preds.multi$aggr$mean - df.test$semimajoraxis)^2))

pdf("pred_parallel.pdf") 
plot(preds.multi$aggr$mean, df.test$semimajoraxis)
dev.off()


round(c(rmse.default, rmse.P50, rmse.parallel),2)


###############################


#Prediction based on the best model () or the MPM (Median Probability Model)

get.best.model(result = result.default)
preds.best <- predict(get.best.model(result.default), df.test[, -1])
sqrt(mean((preds.best - df.test$semimajoraxis)^2))

get.mpm.model(result = result.default, y = df.train$semimajoraxis, x = df.train[, -1])
preds.mpm <- predict(get.mpm.model(result.default, y = df.train$semimajoraxis, x = df.train[, -1]), df.test[, -1])
sqrt(mean((preds.mpm - df.test$semimajoraxis)^2))


####################################################
#
# Diagnostic plots  (Section 3.5)
#
####################################################


pdf("diagn_default.pdf") 
diagn_plot(result.default, ylim = c(600,1500), FUN = max)
dev.off()
diagn_plot(result.default, ylim = c(600,1500), FUN = max)


pdf("diagn_long.pdf") 
diagn_plot(result.P50, ylim = c(600,1500), FUN = max)
dev.off()
diagn_plot(result.P50, ylim = c(600,1500), FUN = max)


pdf("diagn_par.pdf") 
diagn_plot(result_parallel, ylim = c(600,1500),FUN = max)
dev.off()

diagn_plot(result_parallel, ylim = c(600,1500),FUN = max)


