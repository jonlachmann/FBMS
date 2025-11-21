library(microbenchmark)
library(Rcpp)
library(RcppArmadillo)
library(fastglm)
library(BAS)

# Original R function
estimate.logic.tcch.general <- function(y, x, model, complex, params) {
  
  if (length(params) == 0) 
    params <- list(r = 1/dim(x)[1]) 
  
  suppressWarnings({
    mod <- fastglm(as.matrix(x[, model]), y, family = gaussian())
  })
  
  # Compute the general complexity prior
  log_prior <- log(params$r) * sum(complex$oc)
  
  # Compute other terms for tCCH prior on g
  p.v <- (params$n + 1) / (mod$rank + 1)
  
  y_mean <- mean(y)
  TSS <- sum((y - y_mean)^2)
  RSS <- sum(mod$residuals^2)
  R.2 <- 1 - (RSS / TSS)
  p <- mod$rank
  
  # Marginal log-likelihood using tCCH prior for g
  mloglik = (-0.5 * p * log(p.v) 
             - 0.5 * (params$n - 1) * log(1 - (1 - 1/p.v) * R.2) 
             + log(beta((params$p.a + p) / 2, params$p.b / 2)) 
             - log(beta(params$p.a / 2, params$p.b / 2)) 
             + log(phi1(params$p.b / 2, 
                        (params$n - 1) / 2, 
                        (params$p.a + params$p.b + p) / 2, 
                        params$p.s / (2 * p.v), 
                        R.2 / (p.v - (p.v - 1) * R.2))) 
             - hypergeometric1F1(params$p.b / 2, 
                                 (params$p.a + params$p.b) / 2, 
                                 params$p.s / (2 * p.v), log = T)) 
  
  # Stability check
  if (is.na(mloglik) || is.nan(mloglik) || mloglik == -Inf) {
    mloglik = -10000
  }
  
  logpost <- mloglik + log_prior
  
  # Stability check for final log-posterior
  if (logpost == -Inf) {
    logpost = -10000
  }
  
  return(list(crit = logpost, coefs = mod$coefficients))
}


# Rcpp implementation of the function
cppFunction(depends = "RcppArmadillo", code = '
double compute_log_posterior(NumericVector residuals, int p, int n, double r, double p_a, double p_b, double p_s, NumericVector complexity_oc) {
    // Compute R^2
    double RSS = sum(residuals * residuals);
    double TSS = sum((residuals - mean(residuals)) * (residuals - mean(residuals)));
    double R2 = 1 - (RSS / TSS);
    
    // Compute log prior complexity term
    double log_prior_complexity = log(r) * sum(complexity_oc);
    
    // Compute p.v
    double p_v = (n + 1.0) / (p + 1.0);
    
    // Compute marginal log likelihood
    double mloglik = (-0.5 * p * log(p_v) 
                      - 0.5 * (n - 1) * log(1 - (1 - 1 / p_v) * R2) 
                      + R::lbeta((p_a + p) / 2.0, p_b / 2.0) 
                      - R::lbeta(p_a / 2.0, p_b / 2.0) 
                      + log(R::pbeta(R2 / (p_v - (p_v - 1) * R2), p_b / 2.0, (n - 1) / 2.0, 1, 0)) 
                      - R::pbeta(p_s / (2.0 * p_v), p_b / 2.0, (p_a + p_b) / 2.0, 1, 1));
    
    // Stability check
    if (std::isnan(mloglik) || std::isinf(mloglik)) {
        mloglik = -10000;
    }
    
    double logpost = mloglik + log_prior_complexity + n;
    
    if (std::isinf(logpost)) {
        logpost = -10000;
    }
    
    return logpost;
}')



estimate_logic_tcch_rcpp <- function(y, x, model, complex, params) {
  
  if (length(params) == 0) 
    params <- list(r = 1 / nrow(x)) 
  
  # Fit the model using fastglm
  suppressWarnings({
    mod <- fastglm(as.matrix(x[, model]), y, family = gaussian())
  })
  
  # Call the Rcpp function for log-posterior computation
  logpost <- compute_log_posterior(
    residuals = mod$residuals,
    p = mod$rank,
    n = params$n,
    r = params$r,
    p_a = params$p.a,
    p_b = params$p.b,
    p_s = params$p.s,
    complexity_oc = complex$oc
  )
  
  return(list(crit = logpost, coefs = mod$coefficients))
}

# Generate test data
set.seed(42)
n <- 100
p <- 5
X <- matrix(rnorm(n * p), n, p)
beta <- c(1, -1, 0.5, 0, 0)
y <- X %*% beta + rnorm(n)

params <- list(n = n, p.a = 1, p.b = 1, p.s = 0.1, r = 0.01)  # Prior hyperparameters
model <- c(1, 2, 3)  # Assume we select variables 1, 2, 3
complex <- list(oc = c(1, 1, 1))  # Example complexity for selected model

# Convert for Rcpp function
model_vec <- as.integer(model)  # Convert to unsigned integer vector for Rcpp
complex_vec <- as.numeric(complex$oc)

# Run both implementations
result_r <- estimate.logic.tcch.general(y, X, model, complex, params)
result_rcpp <- estimate_logic_tcch_rcpp(y, X, model, complex, params)

# Check if results match
print("Checking results...")
print(all.equal(result_r$crit, result_rcpp$crit, tolerance = 1e-6))
print(all.equal(result_r$coefs, as.numeric(result_rcpp$coefs), tolerance = 1e-6))

# Benchmarking
bench <- microbenchmark(
  R = estimate.logic.tcch.general(y, X, model, complex, params),
  Rcpp = estimate_logic_tcch_rcpp(y, X, model_vec, complex_vec, params$r, params$n, params$p.a, params$p.b, params$p.s),
  times = 100
)

# Print benchmark results
print(bench)

# Calculate speedup
speedup <- median(bench$time[bench$expr == "R"]) / median(bench$time[bench$expr == "Rcpp"])
print(paste("Speedup: ", round(speedup, 2), "x"))



?BAS::bayesglm.fit(x = X, y = y>mean(y),family = binomial(),coefprior = aic.prior)


library(BAS)

set.seed(42)
n <- 100
p <- 5
X <- matrix(rnorm(n * p), n, p)
beta <- c(1, -1, 0.5, 0, 0)
y <- X %*% beta + rnorm(n)

data <- data.frame(y = y>mean(y), x = X)


suppressWarnings({mod <- bas.glm(y ~ 1+x.1, data = data, betaprior = CCH(alpha = 0.5, beta = as.numeric(nrow(data)), s = 0), method = "deterministic", family = binomial(), modelprior = uniform(), n.models = 2, initprobs = 'eplogp', laplace = T)})
mod$logmarg


result <- tryCatch({
result <- .Call(BAS:::C_glm_deterministic,
      y = as.numeric(y>mean(y)), X = cbind(1,X[,1]),
      Roffset = as.numeric(rep(0, length(y))),
      Rweights = as.numeric(rep(1, length(y))),
      Rprobinit = as.numeric(c(1,0.99)),
      Rmodeldim = as.integer(c(0,0)),
      modelprior = uniform(),
      betaprior = CCH(alpha = 0.5, beta = as.numeric(nrow(data)), s = 0),
      family = binomial(), 
      Rcontrol = glm.control(),
      Rlaplace =  as.integer(T))
lm = result$logmarg
print(lm)
}, error = function(e) {
  warning("Error in C call: ", e$message)
  list(logmarg = NA)
})

for (i in 1:1500) {
  result <- .Call(BAS:::C_glm_deterministic,
                  y = as.numeric(y > mean(y)),
                  X = cbind(1, X),
                  Roffset = as.numeric(rep(0, length(y))),
                  Rweights = as.numeric(rep(1, length(y))),
                  Rprobinit = as.numeric(c(1, 0.99)),
                  Rmodeldim = as.integer(0),
                  modelprior = BAS::uniform(),
                  betaprior = betaprior,
                  family = family,
                  Rcontrol = control,
                  Rlaplace = as.integer(FALSE))
  print(result$logmarg)
}

data <- data.frame(y = y, x = X)
mod <- bas.lm(y ~ 1+x.1, data = data, prior =  "ZS-null", alpha = 4, method = "deterministic", modelprior = uniform(), n.models = 2, initprobs = 'eplogp')
mod$logmarg


result <- .Call(BAS:::C_glm_deterministic,
                y = as.numeric(y>mean(y)), X = cbind(1,X[,1]),
                Roffset = as.numeric(rep(0, length(y))),
                Rweights = as.numeric(rep(1, length(y))),
                Rprobinit = as.numeric(c(1,0.99)),
                Rmodeldim = as.integer(c(0,0)),
                modelprior = uniform(),
                betaprior = CCH(alpha = 0.5, beta = as.numeric(nrow(data)), s = 0),
                family = binomial(), 
                Rcontrol = glm.control(),
                Rlaplace =  as.integer(T))

result <- .Call(BAS:::C_deterministic,
                y = y,
                X = cbind(1,X[,1]),
                as.numeric(rep(1, length(y))),
                as.numeric(c(1,0.99)),
                as.integer(c(0,0)),
                incint = as.integer(TRUE),
                alpha = as.numeric(4),
                method = as.integer(0),
                modelprior = uniform(),
                Rpivot = TRUE,
                Rtol = 1e-7
              )
result$logmarg

method.num <- switch(
  prior,
  "g-prior" = 0,
  "hyper-g" = 1,
  "EB-local" = 2,
  "BIC" = 3,
  "ZS-null" = 4,
  "ZS-full" = 5,
  "hyper-g-laplace" = 6,
  "AIC" = 7,
  "EB-global" = 2,
  "hyper-g-n" = 8,
  "JZS" = 9
)