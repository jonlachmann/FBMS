#######################################################
#
# Example 5 (Section 4.3):
#
# Fractional Polynomials: Depths is set to 1
#
# This is the valid version for the JSS Paper
#
#######################################################


library(FBMS)
use.fbms = FALSE  


df = read.csv2(file = "tests/art.csv",sep = ",",dec = ".")[,c(16,1:3,5:8,10:14)]

summary(df)


#number of observations in the data

n = dim(df)[1] 

#number of covariates

p = dim(df)[2] - 1   


set.seed(040590)



mu = 0.1 + p05(df$x1) + df$x1 + pm05(df$x3) + p0pm05(df$x3) + df$x4a + pm1(df$x5) + p0(df$x6) + df$x8 + df$x10


df$y = rnorm(n =n, mean = mu,sd = 1)


transforms <- c("p0","p2","p3","p05","pm05","pm1","pm2","p0p0","p0p05","p0p1","p0p2","p0p3","p0p05","p0pm05","p0pm1","p0pm2")
probs <- gen.probs.gmjmcmc(transforms)
probs$gen <- c(0,1,0,1) # Only modifications!
params <- gen.params.gmjmcmc(df)
params$feat$D <- 1   # Set depth of features to 1

#to set variance to unknown uncomment below
#params$loglik$var <- "unknown"

####################################################
#
# single thread analysis
#
####################################################

set.seed(123)

if (use.fbms) {
  result <- fbms(data = df, method = "gmjmcmc", transforms = transforms, 
                 probs = probs, params = params)
} else {
  result <- gmjmcmc(df, transforms = transforms, probs = probs, params = params)
}
summary(result, labels = names(df[-1]))
#plot.diagn(result,FUN = median)



####################################################
#
# multiple thread analysis
#
####################################################

set.seed(101)

if (use.fbms) {
  result_parallel <- fbms(data = df, method = "gmjmcmc.parallel", transforms = transforms, 
                 probs = probs, params = params, P=25)
} else {
  result_parallel =  gmjmcmc.parallel(runs = 40, cores = 10, data = df, 
                        transforms = transforms, probs = probs, params = params, P=25)
}
#summary(result_parallel, labels = names(df[-1]))
summary(result_parallel, labels = names(df[-1]), tol = 0.01)

diagn_plot(result_parallel, FUN = median)




set.seed(102)

if (use.fbms) {
  result_parallel2 <- fbms(runs = 40, cores = 10,data = df, method = "gmjmcmc.parallel", transforms = transforms, 
                          probs = probs, params = params, P=25, N.init=1000, N.final=2000)
} else {
  result_parallel2 =  gmjmcmc.parallel(runs = 40, cores = 10,data = df, 
                          transforms = transforms, probs = probs, params = params, 
                          P=25, N.init=1000, N.final=2000)
}
#summary(result_parallel2, labels = names(df[-1]))
summary(result_parallel2, labels = names(df[-1]), tol = 0.01)

diagn_plot(result_parallel2,FUN = median)


# Very large number of mjmcmc iterations (not needed for paper)
set.seed(104)

if (use.fbms) {
  result_parallel3 <- fbms(data = df, method = "gmjmcmc.parallel", transforms = transforms, 
                           probs = probs, params = params, P=50, N.init=2000, N.final=4000)
} else {
  result_parallel3 =  gmjmcmc.parallel(runs = 40, cores = 40, data = df, transforms = transforms, 
                                       probs = probs, params = params, P=50, N.init=2000, N.final=4000)
}
#summary(result_parallel3, labels = names(df[-1]))
summary(result_parallel3, labels = names(df[-1]), tol = 0.01)





