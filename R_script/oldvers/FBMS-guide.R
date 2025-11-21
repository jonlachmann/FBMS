ß## ----include=FALSE------------------------------------------------------------
knitr::opts_chunk$set(
  message = TRUE,  # show package startup and other messages
  warning = FALSE, # suppress warnings
  echo    = TRUE,  # show code
  results = "hide" # hide default printed results unless printed via printn()
)

# For careful printing of only what I explicitly ask for
printn <- function(x) {
  # Capture the *exact* console print output as a character vector
  txt <- capture.output(print(x))
  # Combine lines with newline, send as a message to be shown in output
  message(paste(txt, collapse = "\n"))
}

library(FBMS)

## ----eval=FALSE, include=FALSE------------------------------------------------
#  library(FBMS)

## -----------------------------------------------------------------------------

# Parameters for parallel runs are set to a single thread and single core to comply with CRAN requirenments (please tune for your machine if you have more capacity)
runs  <- 1  # 1 set for simplicity; use rather 16 or more
cores <- 1  # 1 set for simplicity; use rather 8 or more

## -----------------------------------------------------------------------------
# Load example
data <- FBMS::exoplanet

# Choose a small but expressive transform set for a quick demo
transforms <- c("sigmoid", "sin_deg", "exp_dbl", "p0", "troot", "p3")

# ---- fbms() call (simple GMJMCMC) ----
# Key parameters (explicit):
# - formula      : semimajoraxis ~ 1 + .       # response and all predictors
# - data         : data                        # dataset
# - beta_prior   : list(type = "g-prior")      # parameter prior    
# - model_prior  : list(r = 1/dim(data)[1])    # model prior   
# - method       : "gmjmcmc"                   # exploration strategy
# - transforms   : transforms                  # nonlinear feature dictionary
# - P            : population size per generation (search breadth)
result_single <- fbms(
  formula     = semimajoraxis ~ 1 + .,
  data        = data,
  beta_prior  = list(type = "g-prior", alpha = dim(data)[1]),
  model_prior = list(r = 1/dim(data)[1]),
  method      = "gmjmcmc",
  transforms  = transforms,
  P           = 20
)

# Summarize
printn(summary(result_single))

## -----------------------------------------------------------------------------

# ---- fbms() call (parallel GMJMCMC) ----
# Key parameters (explicit):
# - formula      : semimajoraxis ~ 1 + .       # response and all predictors
# - data         : data                        # dataset
# - beta_prior   : list(type = "g-prior")      # parameter prior   
# - model_prior  : list(r = 1/dim(data)[1])    # model prior   
# - method       : "gmjmcmc"                   # exploration strategy
# - transforms   : transforms                  # nonlinear feature dictionary
# - runs, cores  : parallelization controls
# - P            : population size per generation (search breadth)
result_parallel <- fbms(
  formula     = semimajoraxis ~ 1 + .,
  data        = data,
  beta_prior  = list(type = "g-prior", alpha = dim(data)[1]),
  model_prior = list(r = 1/dim(data)[1]),
  method      = "gmjmcmc.parallel",
  transforms  = transforms,
  runs        = runs*10,  # by default the rmd has runs =  1; increase for convergence
  cores       = cores,    # by default the rmd has cores = 1; increase for convergence
  P           = 20
)

# Summarize
printn(summary(result_parallel))

## -----------------------------------------------------------------------------
plot(result_parallel)

## -----------------------------------------------------------------------------
diagn_plot(result_parallel)

## -----------------------------------------------------------------------------
library(mvtnorm)

n <- 100               # sample size
p <- 20                # number of covariates
k <- 5                 # size of true submodel

correct.model <- 1:k
beta.k <- (1:5)/5

beta <- rep(0, p)
beta[correct.model] <- beta.k

set.seed(123)
x <- rmvnorm(n, rep(0, p))
y <- x %*% beta + rnorm(n)

# Standardize
y <- scale(y)
X <- scale(x) / sqrt(n)

df <- as.data.frame(cbind(y, X))
colnames(df) <- c("Y", paste0("X", seq_len(ncol(df) - 1)))

printn(correct.model)
printn(beta.k)

## -----------------------------------------------------------------------------
# ---- fbms() call (MJMCMC) ----
# Explicit prior choice:
#   beta_prior = list(type = "g-prior", alpha = 100)
# To switch to another prior, e.g. robust:
#   beta_prior = list(type = "robust")
result.lin <- fbms(
  formula    = Y ~ 1 + .,
  data       = df,
  method     = "mjmcmc",
  N          = 5000,                         # number of iterations
  beta_prior = list(type = "g-prior", alpha = 100)
)

## -----------------------------------------------------------------------------
plot(result.lin)

## -----------------------------------------------------------------------------
# 'effects' specifies quantiles for posterior modes of effects across models
printn(summary(result.lin, effects = c(0.5, 0.025, 0.975)))

## -----------------------------------------------------------------------------
# ---- fbms() call (parallel MJMCMC) ----
# Explicit prior choice:
#   beta_prior = list(type = "g-prior", alpha = 100)
# To switch to another prior, e.g. robust:
#   beta_prior = list(type = "robust")
#   method = mjmcmc.parallel
#   runs, cores  : parallelization controls
result.lin.par <- fbms(
  formula    = Y ~ 1 + .,
  data       = df,
  method     = "mjmcmc.parallel",
  N          = 5000,                         # number of iterations
  beta_prior = list(type = "g-prior", alpha = 100),
  runs = runs,
  cores = cores
)
printn(summary(result.lin.par, effects = c(0.5, 0.025, 0.975)))

## -----------------------------------------------------------------------------
# Create FP-style response with known structure, covariates are from previous example
df$Y <- p05(df$X1) + df$X1 + pm05(df$X3) + p0pm05(df$X3) + df$X4 +
         pm1(df$X5) + p0(df$X6) + df$X8 + df$X10 + rnorm(nrow(df))

# Allow common FP transforms
transforms <- c(
  "p0", "p2", "p3", "p05", "pm05", "pm1", "pm2", "p0p0",
  "p0p05", "p0p1", "p0p2", "p0p3", "p0p05", "p0pm05", "p0pm1", "p0pm2"
)

# Generation probabilities — here only modifications and mutations
probs <- gen.probs.gmjmcmc(transforms)
probs$gen <- c(0, 1, 0, 1)

# Feature-generation parameters
params <- gen.params.gmjmcmc(ncol(df) - 1)
params$feat$D <- 1   # max depth 1 features

## -----------------------------------------------------------------------------
result <- fbms(
  formula    = Y ~ 1 + .,
  data       = df,
  method     = "gmjmcmc",
  transforms = transforms,
  beta_prior = list(type = "Jeffreys-BIC"),
  probs      = probs,
  params     = params,
  P          = 25
)

printn(summary(result))

## -----------------------------------------------------------------------------
result_parallel <- fbms(
  formula    = Y ~ 1 + .,
  data       = df,
  method     = "gmjmcmc.parallel",
  transforms = transforms,
  beta_prior = list(type = "Jeffreys-BIC"),
  probs      = probs,
  params     = params,
  P          = 25,
  runs       = runs,
  cores      = cores
)

printn(summary(result_parallel))

## -----------------------------------------------------------------------------
# Custom approximate log marginal likelihood for mixed model using Laplace approximation
mixed.model.loglik.lme4 <- function (y, x, model, complex, mlpost_params) {
  if (sum(model) > 1) {
    x.model <- x[, model]
    data <- data.frame(y, x = x.model[, -1], dr = mlpost_params$dr)
    mm <- lmer(as.formula(paste0("y ~ 1 +",
                          paste0(names(data)[2:(ncol(data)-1)], collapse = "+"),
                          " + (1 | dr)")), data = data, REML = FALSE)
  } else {
    data <- data.frame(y, dr = mlpost_params$dr)
    mm <- lmer(y ~ 1 + (1 | dr), data = data, REML = FALSE)
  }
  # log marginal likelihood (Laplace approx) + log model prior
  mloglik <- as.numeric(logLik(mm)) - 0.5 * log(length(y)) * (ncol(data) - 2)
  if (length(mlpost_params$r) == 0) mlpost_params$r <- 1 / nrow(x)
  lp <- log_prior(mlpost_params, complex)
  list(crit = mloglik + lp, coefs = fixef(mm))
}

## -----------------------------------------------------------------------------
library(lme4)
data(Zambia, package = "cAIC4")

df <- as.data.frame(sapply(Zambia[1:5], scale))

transforms <- c("p0","p2","p3","p05","pm05","pm1","pm2",
                "p0p0","p0p05","p0p1","p0p2","p0p3",
                "p0p05","p0pm05","p0pm1","p0pm2")

probs <- gen.probs.gmjmcmc(transforms)
probs$gen <- c(1, 1, 0, 1) # include modifications and interactions

params <- gen.params.gmjmcmc(ncol(df) - 1)
params$feat$D <- 1
params$feat$pop.max <- 10

result2a <- fbms(
  formula      = z ~ 1 + .,
  data         = df,
  method       = "gmjmcmc.parallel",
  transforms   = transforms,
  probs        = probs,
  params       = params,
  P            = 25,
  N            = 100,
  runs         = runs,
  cores        = cores,
  family       = "custom",
  loglik.pi    = mixed.model.loglik.lme4,
  model_prior  = list(r = 1 / nrow(df)), # model_prior is passed to mlpost_params
  extra_params = list(dr = droplevels(Zambia$dr)) # extra_params are passed to mlpost_params
)

printn(summary(result2a, tol = 0.05, labels = names(df)[-1]))

## -----------------------------------------------------------------------------
n <- 2000
p <- 50

set.seed(1)
X2 <- as.data.frame(matrix(rbinom(n * p, size = 1, prob = runif(n * p, 0, 1)), n, p))
y2.Mean <- 1 + 7*(X2$V4*X2$V17*X2$V30*X2$V10) + 9*(X2$V7*X2$V20*X2$V12) + 
           3.5*(X2$V9*X2$V2) + 1.5*(X2$V37)

Y2 <- rnorm(n, mean = y2.Mean, sd = 1)
df <- data.frame(Y2, X2)

# Train/test split
df.training <- df[1:(n/2), ]
df.test     <- df[(n/2 + 1):n, ]
df.test$Mean <- y2.Mean[(n/2 + 1):n]

## -----------------------------------------------------------------------------
estimate.logic.lm <- function(y, x, model, complex, mlpost_params) {
  suppressWarnings({
    mod <- fastglm(as.matrix(x[, model]), y, family = gaussian())
  })
  mloglik <- -(mod$aic + (log(length(y)) - 2) * (mod$rank)) / 2
  wj <- complex$width
  lp <- sum(log(factorial(wj))) - sum(wj * log(4 * mlpost_params$p) - log(4))
  logpost <- mloglik + lp
  if (logpost == -Inf) logpost <- -10000
  list(crit = logpost, coefs = mod$coefficients)
}

## -----------------------------------------------------------------------------
set.seed(5001)

# Only "not" operator; "or" is implied by De Morgan via "and" + "not"
transforms <- c("not")
probs <- gen.probs.gmjmcmc(transforms)
probs$gen <- c(1, 1, 0, 1) # no projections

params <- gen.params.gmjmcmc(p)
params$feat$pop.max <- 50
params$feat$L <- 15

result <- fbms(
  formula     = Y2 ~ 1 + .,
  data        = df.training,
  method      = "gmjmcmc",
  transforms  = transforms,
  N           = 500,
  P           = 25,
  family      = "custom",
  loglik.pi   = estimate.logic.lm,
  probs       = probs,
  params      = params,
  model_prior = list(p = p)
)

printn(summary(result))

# Extract models
mpm   <- get.mpm.model(result, y = df.training$Y2, x = df.training[,-1],
                       family = "custom", loglik.pi = estimate.logic.lm, params = list(p = 50))
printn(mpm$coefs)

mpm2  <- get.mpm.model(result, y = df.training$Y2, x = df.training[,-1])
printn(mpm2$coefs)

mbest <- get.best.model(result)
printn(mbest$coefs)

## -----------------------------------------------------------------------------
# Correct link is identity for Gaussian
pred       <- predict(result, x = df.test[,-1], link = function(x) x)
pred_mpm   <- predict(mpm,   x = df.test[,-1], link = function(x) x)
pred_best  <- predict(mbest, x = df.test[,-1], link = function(x) x)

# RMSEs
printn(sqrt(mean((pred$aggr$mean - df.test$Y2)^2)))
printn(sqrt(mean((pred_mpm     - df.test$Y2)^2)))
printn(sqrt(mean((pred_best    - df.test$Y2)^2)))
printn(sqrt(mean((df.test$Mean - df.test$Y2)^2)))

# Errors to the true mean (oracle)
printn(sqrt(mean((pred$aggr$mean - df.test$Mean)^2)))
printn(sqrt(mean((pred_best      - df.test$Mean)^2)))
printn(sqrt(mean((pred_mpm       - df.test$Mean)^2)))

# Quick diagnostic plot
plot(pred$aggr$mean, df.test$Y2,
     xlab = "Predicted (BMA)", ylab = "Observed")
points(pred$aggr$mean, df.test$Mean, col = 2)
points(pred_best,      df.test$Mean, col = 3)
points(pred_mpm,       df.test$Mean, col = 4)

