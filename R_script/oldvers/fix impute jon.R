library(mgcv) 
library(foreign)
#install.packages("XLConnect")
#library(XLConnect)
library(readxl)
# install.packages("e1071")
library(e1071)

#install.packages("knitr")
library(knitr)
library(foreign)
library(readxl)
## tableone package itself
#install.packages("tableone")
#library(tableone)
## plotting
library(ggplot2)
#install.packages("pROC")
library(pROC)
#install.packages("memisc")
#library(memisc)
library(sjlabelled)

#library(gmodels)  # Für CrossTable
library(plyr)
library(dplyr)
#install.packages("devtools")
library(devtools)
library(survival)
## Load rms package
library(rms)

#library(DescTools)

#install.packages("Hmisc")
library(Hmisc)
#install.packages("Gmisc")
#library(Gmisc)
library(corrplot)

library(tidyr)

library(lubridate)
library(mice)
#library("VIM")

library(glmnet)


Excel.File = "/Users/aliaksandrhome/Rprojects/FBMS impute Bones data/Data4.xlsx"

#df <- read.csv2("C:/CC/Projekte/P24110_Kanz/Data3.csv")
#df <- as.data.frame(read_excel(Excel.File))


df <- as.data.frame(read_excel(Excel.File, col_types = c(
  "text", "numeric", "numeric", "numeric", 
  "numeric", "numeric", "numeric", "numeric", 
  "numeric", "numeric", "numeric", "numeric", 
  "numeric", "numeric", "numeric", "numeric", 
  "numeric", "numeric", "numeric", "numeric", 
  "numeric", "numeric", "numeric", "numeric", 
  "numeric", "numeric", "numeric", "numeric", 
  "numeric", "numeric", "numeric", "numeric", 
  "numeric", "numeric", "numeric", "numeric", 
  "numeric", "numeric", "numeric", "numeric", 
  "numeric", "numeric", "numeric", "numeric", 
  "numeric", "numeric", "numeric", "numeric", 
  "numeric", "numeric", "numeric", "numeric", 
  "numeric", "numeric", "numeric", "numeric", 
  "numeric", "numeric", "numeric", "numeric", 
  "numeric", "numeric", "numeric", "numeric", 
  "numeric", "numeric", "numeric", "numeric", 
  "numeric", "numeric", "numeric")))
str(df)

#df$Bef_Nr = as.factor(df$Bef_Nr)
df$TestData = is.na(df$Sex_Molekular)



df.Training = df[df$TestData == F,]
df.Training = df.Training %>% dplyr::select(- C14_mean)

df.Training = df.Training %>% dplyr::select(- c(Clavicula_4, Clavicula_5,Femur_18, Tibia_71,    Tibia_74,Fibula_1,  Fibula_4a ))
m = dim(df.Training)[2]

x.indx = 7:(m-1)

MV = rowSums(!is.na(df.Training[,x.indx[-1]]))

df.Training = df.Training[MV>1,]
n = dim(df.Training)[1]


devtools::load_all()


transforms <- c("p0","p2","p3","p05","pm05","pm1","pm2")
#transforms <- c("p0","p2","p3","p05","pm05","pm1","pm2","p0p0","p0p05","p0p1","p0p2","p0p3","p0p05","p0pm05","p0pm1","p0pm2")

probs <- gen.probs.gmjmcmc(transforms)
probs$gen <- c(1,1,0,1)

df.Training$Bef_Nr <- as.integer(df.Training$Bef_Nr)
df.Training$Sex_Molekular = as.integer(df.Training$Sex_Molekular) - 1

df.test <- df.Training[201:285,]
df.Training <- df.Training[1:200,]

#df.Training$TestData <- NULL

params <- gen.params.gmjmcmc(ncol(df.Training)- 1 + sum(sapply(names(df.Training)[-1],function(name)sum(is.na(df.Training[[name]]))>0)))

params$feat$pop.max <- 150


result_parallel =  fbms(formula = Sex_Molekular ~ 1 + .,runs = 8, cores = 8, data = df.Training, transforms = transforms, probs = probs,params = params, P=25,impute = T, beta_prior = list(type = "Jeffreys-BIC"), method = "gmjmcmc.parallel")

pred.obj = predict(result_parallel, x = df.test[,-2], link = function(x)(1/(1+exp(-x))),x_train = df.Training[,-2])

auc(df.test$Sex_Molekular, pred.obj$aggr$mean)
auc(df.Training$Sex_Molekular, df.Training$Sex_Morpho)
auc(df.test$Sex_Morpho, pred.obj$aggr$mean)


bm <- get.best.model(result_parallel)
pred = predict(object = bm, x = df.test[,-2], link = function(x)(1/(1+exp(-x))),x_train = df.Training[,-2])

auc(df.test$Sex_Molekular, pred)
auc(df.Training$Sex_Molekular, df.Training$Sex_Morpho)
auc(df.test$Sex_Morpho, pred)


mpm <- get.mpm.model(result_parallel,family = "binomial",y = df.Training$Sex_Molekular,x = df.Training[,-2])
pred = predict(object = mpm, x = df.test[,-2], link = function(x)(1/(1+exp(-x))),x_train = df.Training[,-2])

auc(df.test$Sex_Molekular, pred)
auc(df.Training$Sex_Molekular, df.Training$Sex_Morpho)
auc(df.test$Sex_Morpho, pred)




