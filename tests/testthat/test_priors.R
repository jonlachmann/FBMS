# Title     : General tests
# Objective : Test the code
# Created by: jonlachmann
# Created on: 2021-02-25

context("Priors")

test_that("Test various priors through the fbms function", {
  set.seed(123)
  x <- matrix(rnorm(300), 100)
  y <- rnorm(100, 0, 0.5) + rowSums(x[, 1:2])
  expect_true(length(y)==100)
  # validate.model <- function (model, x, y) {
  #   expect_true(all(c(model$marg.probs[1:2] > 0.9, model$marg.probs[3] < 0.9)))
  #   summary <- summary(model, labels = c("a", "b", "c"), tol = -1)
  #   expect_true(all(c(summary$marg.probs[1:2] > 0.9, summary$marg.probs[3] < 0.9)))
  #   plot(model)
  #   pred <- predict(model, x)
  #   # Handle paralell runs
  #   if (!is.null(pred$aggr)) {
  #     pred <- pred$aggr
  #   }
  #   rmse <- sqrt(mean((pred$mean - y)^2))
  #   expect_true(rmse < 0.6)
  #   best_model <- get.best.model(model)
  #   mpm_model <- get.mpm.model(model, y, x)
  # }
  # 
  # # No intercept
  # data <- as.data.frame(cbind(y, x))
  # colnames(data) <- c("y", "a", "b", "c")
  # family <- "gaussian"
  # beta_prior <- list(type = "g-prior", g = 5, a = 3, b = 1, s = 1, rho = 0, v = 1, k = 1)
  # 
  # gaussian_priors <- c(
  #   "CH", "tCCH", "TG", "beta.prime", "intrinsic", "ZS-adapted", "uniform","Jeffreys", "benchmark", "robust",
  #   "g-prior", "hyper-g", "EB-local", "ZS-null", "ZS-full", "BIC", "hyper-g-laplace", "AIC", "EB-global", "hyper-g-n", "JZS",
  #   "Jeffreys-BIC", "g-prior"
  # )
  # 
  # expected <- list(
  #   CH = list(crit = 1.79116859063856e+28, coefs = c(X1 = 0.971122645911043, X2 = 1.02417531317894, X3 = -0.0301201893345434)),
  #   tCCH = list(crit = 1.79116859063856e+28, coefs = c(X1 = 0.971122645911043, X2 = 1.02417531317894, X3 = -0.0301201893345434)),
  #   TG = list(crit = 3.30010826314369e+26, coefs = c(X1 = 0.971122645911043, X2 = 1.02417531317894, X3 = -0.0301201893345434)),
  #   beta.prime = list(crit = 2.37584282903985e+33, coefs = c(X1 = 0.971122645911043, X2 = 1.02417531317894, X3 = -0.0301201893345434)),
  #   intrinsic = list(crit = 16.013997004248, coefs = c(X1 = 0.971122645911043, X2 = 1.02417531317894, X3 = -0.0301201893345434)),
  #   `ZS-adapted` = list(crit = 3.08085895954179e+52, coefs = c(X1 = 0.971122645911043, X2 = 1.02417531317894, X3 = -0.0301201893345434)),
  #   uniform = list(crit = 2.6188594083781e+29, coefs = c(X1 = 0.971122645911043, X2 = 1.02417531317894, X3 = -0.0301201893345434)),
  #   Jeffreys = list(crit = 1.96893489476632e+31, coefs = c(X1 = 0.971122645911043, X2 = 1.02417531317894, X3 = -0.0301201893345434)),
  #   benchmark = list(crit = 1.88194647439014e+31, coefs = c(X1 = 0.971122645911043, X2 = 1.02417531317894, X3 = -0.0301201893345434)),
  #   robust = list(crit = 15.4590106422538, coefs = c(X1 = 0.971122645911043, X2 = 1.02417531317894, X3 = -0.0301201893345434)),
  #   `g-prior` = list(crit = 61.0283622188635, coefs = c(X1 = 0.971122645911043, X2 = 1.02417531317894, X3 = -0.0301201893345434)),
  #   `hyper-g` = list(crit = 88.8999586451438, coefs = c(0.971122645911043, 1.02417531317894, -0.0301201893345434)),
  #   `EB-local` = list(crit = 91.5555793778985, coefs = c(0.971122645911043, 1.02417531317894, -0.0301201893345434)),
  #   `ZS-null` = list(crit = 90.2037895063674, coefs = c(0.971122645911043, 1.02417531317894, -0.0301201893345434)),
  #   `ZS-full` = list(crit = 0, coefs = c(0.971122645911043, 1.02417531317894, -0.0301201893345434)),
  #   BIC = list(crit = -170.850296079658, coefs = c(0.971122645911043, 1.02417531317894, -0.0301201893345434)),
  #   `hyper-g-laplace` = list(crit = 88.8449712472898, coefs = c(0.971122645911043, 1.02417531317894, -0.0301201893345434)),
  #   AIC = list(crit = -166.942540800676, coefs = c(0.971122645911043, 1.02417531317894, -0.0301201893345434)),
  #   `EB-global` = list(crit = 91.5555793778985, coefs = c(0.971122645911043, 1.02417531317894, -0.0301201893345434)),
  #   `hyper-g-n` = list(crit = 90.6511779673767, coefs = c(0.971122645911043, 1.02417531317894, -0.0301201893345434)),
  #   JZS = list(crit = 90.9166314790781, coefs = c(0.971122645911043, 1.02417531317894, -0.0301201893345434)),
  #   `Jeffreys-BIC` = list(crit = -83.4856401007207, coefs = c(X1 = 0.971122645911043, X2 = 1.02417531317894, X3 = -0.0301201893345434))
  # )
  # 
  # for (prior in gaussian_priors) {
  #   if (prior == "TG") next
  #   beta_prior$type <- prior
  #   results[[prior]] <- fbms.mlpost.master(data[, 1], x, c(TRUE, TRUE, TRUE), list(), list(beta_prior = beta_prior, family = "gaussian", r = exp(-0.5)))
  #   print(prior)
  #   print(results[[prior]]$crit)
  #   #expect_equal(results[[prior]]$crit, expected[[prior]]$crit)
  #   #mod1 <- fbms(y ~ ., family = family, beta_prior = beta_prior, method = "mjmcmc", data = data, verbose = FALSE)
  #   #if (!(prior %in% c("hyper-g", "EB-local", "uniform", "ZS-null", "ZS-full", "BIC", "hyper-g-laplace", "AIC", "EB-global", "hyper-g-n", "JZS"))) {
  #   #  validate.model(mod1, x, y)
  #   #}
  # }
  # results2 <- list()
  # for (prior in gaussian_priors) {
  #   results2[[prior]] <- fbms.mlik.master2(data[, 1], x, c(TRUE, TRUE, TRUE), list(), list(prior_beta = prior, family = "gaussian", r = exp(-0.5), g = 5, a = 3, b = 1, s = 1, rho = 0, v = 1, k = 1))
  #   print(prior)
  #   print(results[[prior]]$crit)
  # }
})

test_that("likelihood functions handle rank-deficient design matrices without NA coefficients", {
  set.seed(42)
  n <- 50
  x1 <- rnorm(n)
  # Collinear column (x3 == x1)
  X <- cbind(Intercept = 1, x1 = x1, x2 = rnorm(n), x3 = x1)
  y_gauss <- 1 + 2 * x1 + rnorm(n, sd = 0.5)
  y_binom <- rbinom(n, 1, plogis(0.5 + x1))
  model <- rep(TRUE, ncol(X))
  complex <- list(width = 1, oc = 1)

  res_gauss <- gaussian.loglik(y_gauss, X, model, complex, list(r = 1 / n))
  expect_false(anyNA(res_gauss$coefs))
  expect_true(is.finite(res_gauss$crit))

  res_alpha <- gaussian.loglik.alpha(y_gauss, X, model, complex, list(r = 1 / n))
  expect_false(anyNA(res_alpha$coefs))

  res_logis <- logistic.loglik(y_binom, X, model, complex, list(r = 1 / n, subs = 1))
  expect_false(anyNA(res_logis$coefs))
  expect_true(is.finite(res_logis$crit))
})


