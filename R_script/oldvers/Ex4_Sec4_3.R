#######################################################
#
# Example 4 (Section 4.3):
#
# Fractional Polynomials: Depths is set to 1, using only fbms
#
# This is the valid version for the JSS Paper
#
#######################################################


library(FBMS)


url <- "https://www.uniklinik-freiburg.de/fileadmin/mediapool/08_institute/biometrie-statistik/Dateien/Studium_und_Lehre/Lehrbuecher/Multivariable_Model-building/ART.zip"
temp_dir <- tempfile()
download.file(url, tf <- tempfile(fileext = ".zip"), mode = "wb")
unzip(tf, exdir = temp_dir)

df <- read.csv(file.path(temp_dir, "ART/art", "art.csv"))[,c(16,1:3,5:8,10:14)]

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
params <- gen.params.gmjmcmc(ncol(df) - 1)
params$feat$D <- 1   # Set depth of features to 1


####################################################
#
# single thread analysis
#
####################################################

set.seed(123)
result <- fbms(data = df, method = "gmjmcmc", transforms = transforms, 
                 probs = probs, params = params)
summary(result)



####################################################
#
# multiple thread analysis
#
####################################################

set.seed(101)
result_parallel <- fbms(data = df, method = "gmjmcmc.parallel", transforms = transforms, 
                          probs = probs, params = params, P=25,runs = 40, cores = 40)
summary(result_parallel, tol = 0.05)

diagn_plot(result_parallel, FUN = median)




set.seed(102)
  result_parallel2 <- fbms(data = df, method = "gmjmcmc.parallel", transforms = transforms, 
                           probs = probs, params = params, P=25, N=1000, N.final=2000, 
                           runs = 40, cores = 40,)

summary(result_parallel2, tol = 0.05)

diagn_plot(result_parallel2,FUN = median)


##########################

# Using Jeffreys-BIC prior

set.seed(103)
result_parallel3 <- fbms(data = df, method = "gmjmcmc.parallel", beta_prior = list(type = "Jeffreys-BIC"), transforms = transforms, 
                         probs = probs, params = params, P=25, N=1000, N.final=2000, 
                         runs = 40, cores = 40,)

summary(result_parallel3, tol = 0.05)

