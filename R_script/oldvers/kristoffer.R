
library(BAS)

x <- read.csv("https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/refs/heads/master/supplementaries/Mode%20Jumping%20MCMC/supplementary/examples/US%20Data/simcen-x1.txt",header = F)
y <- read.csv(
  "https://raw.githubusercontent.com/aliaksah/EMJMCMC2016/refs/heads/master/supplementaries/Mode%20Jumping%20MCMC/supplementary/examples/US%20Data/simcen-y1.txt",header = F)



data <- data.frame(x,y)

names(data)[16] = "Crime"


res <- BAS::bas.lm(formula = Crime ~.,data = data,prior = "g-prior",alpha = 47,n.models = 32768,method = "deterministic", modelprior = uniform())

strue <- summary(res)

strue[1+c(8,13,14,12,5,9,7,4,6,1,3,2,11,10,15),1]

res.bas1 <- BAS::bas.lm(formula = Crime ~.,data = data,prior = "g-prior",method = "BAS",alpha = 47,n.models=3276,update=500,modelprior= uniform(),initprobs="Uniform")

sbas1 <- summary(res.bas1)

sbas1[1+c(8,13,14,12,5,9,7,4,6,1,3,2,11,10,15),1]


res.bas2 <- BAS::bas.lm(formula = Crime ~.,data = data,prior = "g-prior",method = "BAS",alpha = 47,n.models=3276,update=NULL,modelprior= uniform(),initprobs="Uniform")

sbas2 <- summary(res.bas2)

sbas2[1+c(8,13,14,12,5,9,7,4,6,1,3,2,11,10,15),1]


res.bas3 <- BAS::bas.lm(formula = Crime ~.,data = data,prior = "g-prior",method = "BAS",alpha = 47,n.models=3276,update=500,modelprior= uniform(),initprobs="eplogp")

sbas3 <- summary(res.bas3)

sbas3[1+c(8,13,14,12,5,9,7,4,6,1,3,2,11,10,15),1]


res.mc31 <- BAS::bas.lm(formula = Crime ~.,data = data,prior = "g-prior",alpha = 47,method = "MCMC",n.models = 327, update=NULL,modelprior= uniform(),initprobs="Uniform")

smc31 <- summary(res.mc31)

smc31[1+c(8,13,14,12,5,9,7,4,6,1,3,2,11,10,15),1]




params <- gen.params.mjmcmc(15)
probs <- gen.probs.mjmcmc()



res100 <- mjmcmc(x = x,y = y,loglik.pi = fbms.mlik.master.temp,mlpost_params =  list(family = "gaussian", beta_prior = list(type = "g-prior", g = 47), temp = 100, r = 0.5),N = 100000,params = params,probs = probs)

res1 <- mjmcmc(x = x, y = y,loglik.pi = fbms.mlik.master.temp,mlpost_params =  list(family = "gaussian", beta_prior = list(type = "g-prior", g = 47), temp = 1, r = 0.5), N = 100000,params = params,probs = probs)

res0.1 <- mjmcmc(x = x,y = y,loglik.pi = fbms.mlik.master.temp,mlpost_params =  list(family = "gaussian", beta_prior = list(type = "g-prior", g = 47), temp = 0.5, r = 0.5), N = 100000,params = params,probs = probs)


results <- data.frame(true = strue[2:16,1], t100 = res100$marg.probs[1,],t1 = res1$marg.probs[1,], t01 = res0.1$marg.probs[1,], BAS1 = sbas1[2:16,1],BAS2 = sbas2[2:16,1],BAS3 = sbas3[2:16,1],mc31 = smc31[2:16,1])[c(8,13,14,12,5,9,7,4,6,1,3,2,11,10,15),]

abs(results$true - results)*100















