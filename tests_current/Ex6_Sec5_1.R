#######################################################
#
# Example 6 (Section 5.1): Sanger data again
#
# High dimensional analysis without nonlinearities, using only FBMS
#
# This is the valid version for the JSS Paper
#
#######################################################

#library(devtools)
#devtools::install_github("jonlachmann/FBMS@v1_arxiv", force=T, build_vignettes=F)

library(FBMS)
library(xtable)
library(tictoc)
run.parallel <- TRUE  # Flag to control whether to run gmjmcmc in parallel or just load results


data(SangerData2)
df <- SangerData2
# Rename columns for clarity: response is "y", predictors "x1", "x2", ..., "xp"
colnames(df) = c("y",paste0("x",1:(ncol(df)-1)))
n = dim(df)[1]; p=dim(df)[2]-1

#Use only linear terms and mutations
transforms = c("")
probs = gen.probs.gmjmcmc(transforms)
probs$gen = c(0,0,0,1)


# Select candidate features for the first MJMCMC round by correlation with response
c.vec = unlist(mclapply(2:ncol(df), function(x)abs(cor(df[,1],df[,x]))))
ids = sort(order(c.vec,decreasing=TRUE)[1:50])

# Generate default parameters for GMJMCMC for p-1 predictors
params = gen.params.gmjmcmc(p)
# Restrict feature pre-filtering to top 50 predictors selected by correlation
params$feat$prel.filter <- ids

params$feat$pop.max <- 50    # Maximum population size for the GMJMCMC search

####################################################
#
# Three independent runs of gmjmcmc.parallel
#
####################################################


if (run.parallel) {
  set.seed(12)
    result_parallel1 = fbms(data=df, transforms=transforms,
                            beta_prior = list(type="g-prior", g=max(n,p^2)),
                            probs=probs,params=params,
                            method="gmjmcmc.parallel",
                            P=50,N=1000,runs=10,cores=10)
  save(result_parallel1,file="Ex6_parallel1_orig.RData")

  set.seed(1234)
     result_parallel2=fbms(data=df, transforms=transforms,
                           beta_prior = list(type="g-prior", g=max(n,p^2)),
                           probs=probs,params=params,
                           method="gmjmcmc.parallel",
                           P=50,N=1000,runs=10,cores=10)
  save(result_parallel2,file="Ex6_parallel2_orig.RData")

  set.seed(123456)
     result_parallel3=fbms(data=df, transforms=transforms,
                           beta_prior = list(type="g-prior", g=max(n,p^2)),
                           probs=probs,params=params,
                           method="gmjmcmc.parallel",
                           P=50,N=1000,runs=10,cores=10)
  save(result_parallel3,file="Ex6_parallel3_orig.RData")

} else {

  # If not running gmjmcmc.parallel again, load previously saved results
  load("Ex6_parallel1_orig.RData")
  load("Ex6_parallel2_orig.RData")
  load("Ex6_parallel3_orig.RData")

}


# Summarize results from each of the three parallel runs with tolerance of 0.01

res1 = summary(result_parallel1,tol=0.05)
res1$marg.probs = round(res1$marg.probs,3)
res2 = summary(result_parallel2,tol=0.05)
res2$marg.probs = round(res2$marg.probs,3)
res3 = summary(result_parallel3,tol=0.05)
res3$marg.probs = round(res3$marg.probs,3)

# Combine unique feature names found in all three runs
names.best = unique(c(res1$feats.strings,res2$feats.strings,res3$feats.strings))

# Find maximum number of rows across summaries to equalize sizes for cbind
m = max(nrow(res1),nrow(res2),nrow(res3))
# Pad shorter summaries with empty rows to make them all length m
while(nrow(res1)<m)
  res1 = rbind(res1,c("",""))
while(nrow(res2)<m)
  res2 = rbind(res2,c("",""))
while(nrow(res3)<m)
  res3 = rbind(res3,c("",""))

# Create Latex Table
names(res1) = c("feats","prob")
names(res2) = c("feats","prob")
names(res3) = c("feats","prob")
foo = cbind(res1[1:m,],res2[1:m,],res3[1:m,])
show(print(xtable(foo),include.rownames=FALSE))

# Plot correlation matrix of selected features to visually inspect multicollinearity
ind = order(as.numeric(substring(names.best,first=2)))
names.best = names.best[ind]
X.best = df[,names.best]
pdf("crossplot_Sanger_best.pdf")
corrplot::corrplot(cor(X.best))
dev.off()

